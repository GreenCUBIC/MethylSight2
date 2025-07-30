import os, sys
import argparse

import Bio.SeqIO as SeqIO
from loguru import logger
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import dotenv
import tqdm
import h5py
from sklearn.metrics import average_precision_score, roc_auc_score

import math
import torch
import torchvision.ops.focal_loss as focal_loss
from torch.utils.data import Dataset, random_split
import torch.nn as nn
import torch.nn.functional as F
import lightning as L
from torch.utils.data.sampler import Sampler


from transformers import T5Tokenizer, T5EncoderModel                                                       

import torchvision.ops.focal_loss as focal_loss

MAX_LOGIT = 5.23

def sigmoid(x):
    return 1/(1 + np.exp(-x))

def embed_sequences(sequences, gpu=True):
    
    tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_bfd', do_lower_case=False) 
    gpu_model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_bfd")
    cpu_model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_bfd")
    
    
    gpu_model.eval()
    cpu_model.eval()
    
    sequence_as_aas = [' '.join(list(seq[1])) for seq in sequences]
    
    ids = tokenizer.batch_encode_plus(sequence_as_aas, add_special_tokens=True, padding=True)
    input_ids = torch.tensor(ids['input_ids'])
    attention_mask = torch.tensor(ids['attention_mask'])
    
    if gpu:
        gpu = torch.device('cuda:0')
        gpu_model.to(gpu)
        try:
            input_ids = input_ids.to(gpu)
            attention_mask = attention_mask.to(gpu)
            with torch.no_grad():
                embeddings = gpu_model(input_ids, attention_mask=attention_mask)
        except torch.OutOfMemoryError:
            input_ids = input_ids.to('cpu')
            attention_mask = attention_mask.to('cpu')
            with torch.no_grad():
                embeddings = cpu_model(input_ids, attention_mask=attention_mask)
    else:
        input_ids = input_ids.to('cpu')
        attention_mask = attention_mask.to('cpu')
        with torch.no_grad():
            embeddings = cpu_model(input_ids, attention_mask=attention_mask)   
            
    embeddings = embeddings.last_hidden_state.cpu()
            
    reps = { }
    for i, sequence in enumerate(sequences):
        reps[sequence[0]] = {}
        sequence_embedding = embeddings[i][:len(sequence[1])]
        site_representations = extract_site_representations(sequence[1], sequence_embedding)
        
        for position, site_representation in site_representations.items():
            reps[sequence[0]][position] = site_representation
            
    return reps
            

MODIFICATIONS = {
    'methylation': 0,
    'acetylation': 1,
    'ubiquitination': 2,
    'sumoylation': 3
}

def pr_at_re(x_hat, x_true, recall=0.5):
    
    df = pd.DataFrame({ 'score': x_hat, 'label': x_true })
    
    thresholds = []
    recalls = [] 
    
    index_past_threshold = -1
    
    for i, threshold in enumerate(np.linspace(df.score.min(), df.score.max(), num=1000)):
        thresholds.append(threshold)
        tp = len(df[(df['score'] >= threshold) & (df['label'] == 1)])
        fp = len(df[(df['score'] >= threshold) & (df['label'] == 0)])
        fn = len(df[(df['score'] < threshold) & (df['label'] == 1)])
        re = tp / (tp + fn)
        recalls.append(re)

        if re < recall:
            index_past_threshold = i
            break
        
    if index_past_threshold == -1:
        return 0
        
    t = thresholds[index_past_threshold - 1] + ((recall - recalls[index_past_threshold - 1]) * (thresholds[index_past_threshold] - thresholds[index_past_threshold - 1]) / (recalls[index_past_threshold] - recalls[index_past_threshold - 1]))
    
    
    tp = len(df[(df['score'] >= threshold) & (df['label'] == 1)])
    fp = len(df[(df['score'] >= threshold) & (df['label'] == 0)])
    fn = len(df[(df['score'] < threshold) & (df['label'] == 1)])
    re = tp / (tp + fn)
    pr = tp / (tp + fp)
    
    return pr
    

class ClassificationHead(nn.Module):
    
    def __init__(self, d_model, window_size, layer_widths, dropout=0.15):
        super().__init__()
        
        layers = []
        input_dims = int(d_model * window_size)
        for i in range(len(layer_widths)):
            layers.append(nn.Sequential(
                nn.Linear(input_dims, layer_widths[i]),
                nn.ReLU(),
                nn.Dropout(dropout),
            ))
            input_dims = layer_widths[i]
            
        layers.append(nn.Sequential(
            nn.Linear(input_dims, 1),
            nn.ReLU(),
        ))
        
        self.layers = nn.Sequential(*layers)

    def forward(self, x):
        return self.layers(x)

class MultitaskSampler(Sampler):
    def __init__(self, data, batch_size) -> None:
        self.data = data.reset_index(drop=True)
        self.batch_size = batch_size
        
        self.methylation_indices = np.array(self.data[self.data['modification'] == 'methylation'].index, dtype=int)
        self.acetylation_indices = np.array(self.data[self.data['modification'] == 'acetylation'].index, dtype=int)
        self.ubiquitination_indices = np.array(self.data[self.data['modification'] == 'ubiquitination'].index, dtype=int)
        self.sumoylation_indices = np.array(self.data[self.data['modification'] == 'sumoylation'].index, dtype=int)
        
        self.num_methylation_batches = (len(self.methylation_indices) + self.batch_size - 1) // self.batch_size
        self.num_acetylation_batches = (len(self.acetylation_indices) + self.batch_size - 1) // self.batch_size
        self.num_ubiquitination_batches = (len(self.ubiquitination_indices) + self.batch_size - 1) // self.batch_size
        self.num_sumoylation_batches = (len(self.sumoylation_indices) + self.batch_size - 1) // self.batch_size
        
    def __len__(self) -> int:
        # number of batches to be sampled
        return self.num_methylation_batches + self.num_acetylation_batches + self.num_ubiquitination_batches + self.num_sumoylation_batches

    def __iter__(self):
        # Group into batches where all instances are of the same task 
        # and yield the batches (steps)
        
        methylation_indices = np.copy(self.methylation_indices)
        acetylation_indices = np.copy(self.acetylation_indices)
        ubiquitination_indices = np.copy(self.ubiquitination_indices)
        sumoylation_indices = np.copy(self.sumoylation_indices)
        
        np.random.shuffle(methylation_indices)
        np.random.shuffle(acetylation_indices)
        np.random.shuffle(ubiquitination_indices)
        np.random.shuffle(sumoylation_indices)
        
        methylation_batches = torch.chunk(torch.IntTensor(methylation_indices), self.num_methylation_batches)
        acetylation_batches = torch.chunk(torch.IntTensor(acetylation_indices), self.num_acetylation_batches)
        ubiquitination_batches = torch.chunk(torch.IntTensor(ubiquitination_indices), self.num_ubiquitination_batches)
        sumoylation_batches = torch.chunk(torch.IntTensor(sumoylation_indices), self.num_sumoylation_batches)
        
        for batch in methylation_batches + acetylation_batches + ubiquitination_batches + sumoylation_batches:
            yield batch.tolist()

class KmeDataset(Dataset):
    def __init__(self, embeddings, dataset, window_size):
        self.sites = dataset
        self.embeddings = embeddings
        self.labels = list(self.sites['label'])
        self.window_size = window_size
        
    def __len__(self):
        return len(self.labels)
    
    def __getitem__(self, idx):
        # Use zero padding (i.e. embedding with all zeros) when the chain is too close to the end of the chain
        
        protein = self.sites.iloc[idx]['protein']
        position = self.sites.iloc[idx]['uniprot_position'] 
        modification = self.sites.iloc[idx]['modification'] 
        position_index = position - 1
        representation = torch.Tensor(np.array(self.embeddings[protein]))
        
        protein_len = representation.shape[0]
        representation_dim = representation.shape[1]
        half_window = int((self.window_size - 1)/2)
        padding_left_required = -1 * min(0, position_index - half_window)
        padding_right_required = max(0, position_index + half_window - protein_len + 1)
        
        if padding_left_required > 0:
            representation = torch.cat([torch.zeros((padding_left_required, representation_dim)), representation], dim=0)
        if padding_right_required > 0:
            representation = torch.cat([representation, torch.zeros((padding_right_required, representation_dim))], dim=0)
            
        representation = representation[position_index + padding_left_required - half_window: position_index + padding_left_required + half_window + 1]
        
        # Prepend task token (IN THIS IMPLEMENTATION OF MULTITASK LEARNING, WE DO NOT DO THIS)
        # representation = torch.cat([TOKEN_VALUE[modification] * torch.ones(1, representation_dim), representation])
        
        label = float(self.labels[idx])
        
        return modification, representation, label
    
    
def extract_site_representations(sequence, sequence_embedding, window_size=31):
    
    lysine_indices = [i for i, aa in enumerate(sequence) if aa == 'K']
    
    representations = {}
    
    for position_index in lysine_indices:
        protein_len = sequence_embedding.shape[0]
        representation_dim = sequence_embedding.shape[1]
        half_window = int((window_size - 1)/2)
        padding_left_required = -1 * min(0, position_index - half_window)
        padding_right_required = max(0, position_index + half_window - protein_len + 1)
        
        representation = sequence_embedding.clone().detach()
        
        if padding_left_required > 0:
            representation = torch.cat([torch.zeros((padding_left_required, representation_dim)), representation], dim=0)
        if padding_right_required > 0:
            representation = torch.cat([representation, torch.zeros((padding_right_required, representation_dim))], dim=0)
            
        representation = representation[position_index + padding_left_required - half_window: position_index + padding_left_required + half_window + 1]
        
        if position_index == 47:
            print(representation)

        representations[position_index + 1] = representation
        
    return representations

class PositionalEncoding(nn.Module):

    def __init__(self, d_model: int, dropout: float = 0.1, max_len: int = 5000):
        super().__init__()
        
        self.dropout = nn.Dropout(p=dropout)

        position = torch.arange(max_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model)) # do i need to modify this?
        pe = torch.zeros(1, max_len, d_model)
        pe[0, :, 0::2] = torch.sin(position * div_term) # modified axes here, was pe[:, 0, 0::2]
        pe[0, :, 1::2] = torch.cos(position * div_term) # modified axes here, was pe[:, 0, 1::2]
        self.register_buffer('pe', pe)
        
    def forward(self, x):
        x = x + self.pe[:x.size(0)]
        return self.dropout(x)

#@torch.compile
class Model(L.LightningModule):
    def __init__(self, hparams):
        super().__init__()
        self.save_hyperparameters()
        self.lr_step_size = hparams['lr_step_size']
        self.lr_gamma = hparams['lr_gamma']
        self.learning_rate = hparams['learning_rate']
        self.batch_size = hparams['batch_size']
        self.methylation_loss_factor = hparams['methylation_loss_factor']
        self.loss_function = eval(hparams['loss_function'])
        self.training_step_outputs = []
        self.validation_step_outputs = []
        
        self.embedding_layer = nn.Sequential(
            nn.Linear(hparams['input_dims'], hparams['embedder_width']),
            nn.Dropout(hparams['dropout']),
            nn.ReLU(),
            nn.Linear(hparams['embedder_width'], hparams['d_model']),
            nn.Dropout(hparams['dropout']),
            nn.ReLU(),
        )

        self.positional_encoder = PositionalEncoding(hparams['d_model'], dropout=hparams['dropout'], max_len=hparams['window_size'])
        
        transformer_layers = []

        for i in range(len(hparams['n_heads'])):
            transformer_layers.append(
                nn.TransformerEncoderLayer(
                    hparams["d_model"],
                    hparams["n_heads"][i],
                    dropout=hparams["dropout"],
                    activation=nn.ReLU()
                ),
            )

        self.transformer_layers = nn.Sequential(*transformer_layers)
        
        self.flatten = nn.Flatten()
        
        self.methylation_head = ClassificationHead(hparams['d_model'], hparams['window_size'], hparams['hidden_layer_widths'], dropout=hparams['dropout'])
        self.acetylation_head = ClassificationHead(hparams['d_model'], hparams['window_size'], hparams['hidden_layer_widths'], dropout=hparams['dropout'])
        self.ubiquitination_head = ClassificationHead(hparams['d_model'], hparams['window_size'], hparams['hidden_layer_widths'], dropout=hparams['dropout'])
        self.sumoylation_head = ClassificationHead(hparams['d_model'], hparams['window_size'], hparams['hidden_layer_widths'], dropout=hparams['dropout'])
        
    def forward(self, x, task):
        x = self.embedding_layer(x)
        x = self.positional_encoder(x)
        x = self.transformer_layers(x) 
        x = self.flatten(x)
        
        if task == 'methylation':
            logits = self.methylation_head(x)
        elif task == 'acetylation':
            logits = self.acetylation_head(x)
        elif task == 'ubiquitination':
            logits = self.ubiquitination_head(x)
        elif task == 'sumoylation':
            logits = self.sumoylation_head(x)
        else:
            raise f"Invalid task `{task}` provided."
            
        return logits
    
    def training_step(self, batch, batch_idx):
        tasks, x, y = batch
        task = tasks[0]
        logits = self.forward(x, task)
        y_cpu = y.cpu().detach().numpy()
        y_hat = logits.squeeze(-1)
        y_hat_cpu = y_hat.cpu().detach().numpy()
        loss = self.loss_function(y_hat, y).cpu()
        
        if task == 'methylation':
            loss *= self.methylation_loss_factor
        
        self.log('metrics/batch/loss', loss)
        metrics = { 'loss': loss, 'y': y_cpu, 'y_hat': y_hat_cpu }
        self.training_step_outputs.append(metrics)
        return metrics
    
    def on_training_epoch_end(self):
        loss = np.array([])
        y = np.array([])
        y_hat = np.array([])
        for results_dict in self.training_step_outputs:
            loss = np.append(loss, results_dict["loss"])
            y = np.append(y, results_dict["y"])
            y_hat = np.append(y_hat, results_dict["y_hat"])
        auprc = average_precision_score(y, y_hat)
        self.log("metrics/epoch/loss", loss.mean())
        self.log("metrics/epoch/auprc", auprc)
        self.training_step_outputs.clear()

    def validation_step(self, batch, batch_idx):
        task, x, y = batch
        logits = self.forward(x, task[0])
        y_cpu = y.cpu().detach().numpy()
        y_hat = logits.squeeze(-1)
        y_hat_cpu = y_hat.cpu().detach().numpy()
        loss = self.loss_function(y_hat, y).cpu()
        metrics = { 'loss': loss, 'y': y_cpu, 'y_hat': y_hat_cpu }
        self.validation_step_outputs.append(metrics)
        
    def on_validation_epoch_end(self):
        loss = np.array([])
        y = np.array([])
        y_hat = np.array([])
        for results_dict in self.validation_step_outputs:
            loss = np.append(loss, results_dict["loss"])
            y = np.append(y, results_dict["y"])
            y_hat = np.append(y_hat, results_dict["y_hat"])
        auprc = average_precision_score(y, y_hat)
        auroc = roc_auc_score(y, y_hat)
        prat50re = pr_at_re(y_hat, y, recall=0.5)
        self.logger.experiment["val/loss"] = loss.mean()
        self.logger.experiment["val/auprc"] = auprc
        self.logger.experiment["val/auroc"] = auroc
        self.logger.experiment["val/prat50re"] = prat50re
        self.log("val/loss", loss.mean())
        self.log("val/auprc", auprc)
        self.log("val/auroc", auroc)
        self.log("val/prat50re", prat50re)
        self.validation_step_outputs.clear()
        
    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.learning_rate)
        lr_scheduler = {
            'scheduler': torch.optim.lr_scheduler.StepLR(optimizer, step_size=self.lr_step_size, gamma=self.lr_gamma),
            'name': 'linear_scheduler'
        }
       
        return [optimizer], [lr_scheduler]
 
 
if __name__ == '__main__': 
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, type=str, help='Path to a FASTA file with sequences.')
    parser.add_argument('-w', '--weights', required=True, type=str, help='Path to model checkpoints.')
    parser.add_argument('-o', '--output', required=True, type=str, help='Path to output file.')
    args = vars(parser.parse_args())

    # Load sequences
    sequences = [(x.id, str(x.seq)) for x in SeqIO.parse(args['input'], 'fasta')]

    # Load model
    model = Model.load_from_checkpoint(args['weights'])
    model.eval()

    # Perform inference
    if torch.cuda.is_available():
        model.to('cuda')
        
        scores = []
        logger.info("Embedding the sequences on the GPU...")
        site_embeddings = embed_sequences(sequences, gpu=True)
        logger.info("Making predictions on the GPU...")
        
        for protein, sites in tqdm.tqdm(site_embeddings.items()):
            try:
                for position, representation in sites.items():
                    emb = representation.to('cuda')
                    with torch.no_grad():
                        logit = model(emb, 'methylation').cpu().squeeze(-2).detach().numpy()[0]
                        scores.append({
                            'protein': protein,
                            'position': position,
                            #'logit': float(logit),
                            #'uncorrected_score': float(sigmoid(logit)),
                            'score': float(sigmoid(logit - MAX_LOGIT))
                        })
            except Exception as e:
                print(e)
                print(f"Could not do {protein}... skipping.")
                continue
                    
        df = pd.DataFrame(scores)
        df.to_csv(args['output'], index=False)

    else:

        scores = []
        logger.info("Embedding the sequences on the CPU...")
        site_embeddings = embed_sequences(sequences, gpu=False)
        logger.info("Making predictions on the CPU...")

        for protein, sites in tqdm.tqdm(site_embeddings.items()):
            try:
                for position, representation in sites.items():
                    emb = representation
                    with torch.no_grad():
                        logit = model(emb, 'methylation').squeeze(-2)[0]
                        scores.append({
                            'protein': protein,
                            'position': position,
                            #'logit': float(logit),
                            #'uncorrected_score': float(sigmoid(logit)),
                            'score': float(sigmoid(logit - MAX_LOGIT))
                        })
            except Exception as e:
                print(e)
                print(f"Could not do {protein}... skipping.")
                continue

        df = pd.DataFrame(scores)
        df.to_csv(args['output'], index=False)
