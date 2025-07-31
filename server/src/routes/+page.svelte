<script>
  import { Spinner, Alert, AccordionItem, Accordion } from "flowbite-svelte";
  import {
    Table,
    TableBody,
    TableBodyCell,
    TableBodyRow,
    TableHead,
    TableHeadCell,
    Checkbox,
    TableSearch,
  } from "flowbite-svelte";
  import { Plot, Line, Dot, AxisX, AxisY } from "svelteplot";
  import ms2performance from "../methylsight2_performance_1to36.json";
  import Slider from "@bulatdashiev/svelte-slider";

  let CONSERVATIVE_THRESHOLD = 0.897;
  let PERMISSIVE_THRESHOLD = 0.77;

  let proteinName = $state(null);
  let results = $state(null);

  //let results = null;
  let sequenceText = $state("");
  let success = $state(false);
  let accessionIds = $state("");
  let file = $state(null);
  let isSubmitting = $state(false);
  let message = $state("");
  let threshold = $state([0.9]);
  let performanceAtThreshold = $derived(
    ms2performance[`${threshold[0].toFixed(3)}`]
  );

  const handleFileChange = (event) => {
    file = event.target.files[0];
  };

  const handleSubmit = async (event) => {
    event.preventDefault();
    isSubmitting = true;
    message = "";

    const formData = new FormData();
    formData.append("sequenceText", sequenceText);
    formData.append("accessionIds", accessionIds);
    if (file) {
      formData.append("file", file);
    }

    try {
      const res = await fetch("/upload", {
        method: "POST",
        body: formData,
      });

      const resp = await res.json();
      proteinName = resp.proteinName;
      results = resp.predictions;
      message = resp.message;
      success = true;
    } catch (error) {
      console.error(error);
      message = "An error occurred while uploading.";
    } finally {
      isSubmitting = false;
    }
  };

  const downloadPlainTextCitation = (type) => {
    let data;
    if (type === "bib") {
      data = `@misc{Charih2025,
  title = {Sequence-Based Protein-Protein Interaction Prediction and Its Applications in Drug Discovery},
  author = {Charih, Fran{\\c c}ois and Green, James R. and Biggar, Kyle K.},
  year = {2025},
  publisher = {arXiv},
  doi = {10.48550/ARXIV.2507.19805},
  abstract = {Aberrant protein-protein interactions (PPIs) underpin a plethora of human diseases, and disruption of these harmful interactions constitute a compelling treatment avenue. Advances in computational approaches to PPI prediction have closely followed progress in deep learning and natural language processing. In this review, we outline the state-of the-art for sequence-based PPI prediction methods and explore their impact on target identification and drug discovery. We begin with an overview of commonly used training data sources and techniques used to curate these data to enhance the quality of the training set. Subsequently, we survey various PPI predictor types, including traditional similarity-based approaches, and deep learning-based approaches with a particular emphasis on the transformer architecture. Finally, we provide examples of PPI prediction in systems-level proteomics analyses, target identification, and design of therapeutic peptides and antibodies. We also take the opportunity to showcase the potential of PPI-aware drug discovery models in accelerating therapeutic development.},
  copyright = {arXiv.org perpetual, non-exclusive license},
  keywords = {Biomolecules (q-bio.BM),FOS: Biological sciences,FOS: Computer and information sciences,Machine Learning (cs.LG)},
}`;
    } else {
      data =
        'Charih, François, James R. Green, and Kyle K. Biggar. "Sequence-based Protein-protein Interaction Prediction and Its Applications in Drug Discovery." ArXiv, (2025). Accessed July 31, 2025. https://arxiv.org/abs/2507.19805.';
    }

    const blob = new Blob([data], { type: "text/plain" });
    const fileURL = URL.createObjectURL(blob);
    const downloadLink = document.createElement("a");
    downloadLink.href = fileURL;
    downloadLink.download = "Charih2025." + type;
    document.body.appendChild(downloadLink);
    downloadLink.click();
    URL.revokeObjectURL(fileURL);
  };

  const downloadAsCSV = () => {
    let csv =
      "Protein,Position,MethylSight2 score,Methylated (at Pr=0.75),Methylated (at Re=0.50)";

    results.forEach((r) => {
      csv =
        csv +
        `\n${proteinName},${r.position},${r.score},${r.score > CONSERVATIVE_THRESHOLD ? "Yes" : "No"},${r.score > PERMISSIVE_THRESHOLD ? "Yes" : "No"}`;
    });

    const blob = new Blob([csv], { type: "application/csv" });
    const fileURL = URL.createObjectURL(blob);
    const downloadLink = document.createElement("a");
    downloadLink.href = fileURL;
    downloadLink.download = `MS2_predictions.csv`;
    document.body.appendChild(downloadLink);
    downloadLink.click();
    URL.revokeObjectURL(fileURL);
  };
</script>

<main class="container">
  <section class="hero">
    <h1>MethylSight 2.0</h1>
    <h2>François Charih, Mullen Boulter, Kyle K. Biggar, James R. Green</h2>
    <p>A fully sequence-based lysine methylation predictor</p>
  </section>

  <Accordion flush>
    <AccordionItem open={!results}>
      {#snippet header()}Sequence submission{/snippet}
      <div>
        <form on:submit|preventDefault={handleSubmit} class="upload-form">
          <label for="accessionIds"
            >Provide UniProt accession IDs <b
              >(human proteins in Swiss-Prot only; separated by commas)</b
            >:</label
          >
          <textarea
            id="accessionIds"
            bind:value={accessionIds}
            rows="1"
            placeholder="Q9H7B4"
          ></textarea>

          <label for="sequenceText"
            ><b>or</b> enter a valid protein sequence (FASTA format):</label
          >
          <textarea
            id="sequenceText"
            bind:value={sequenceText}
            rows="8"
            placeholder=">sp|P12345|PROT_HUMAN...\nMSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQ..."
          ></textarea>

          <label for="file"
            ><b>or</b> upload a sequence file (in FASTA format):</label
          >
          <input
            type="file"
            id="file"
            accept=".fasta,.txt"
            on:change={handleFileChange}
          />

          <button type="submit" disabled={isSubmitting}>
            {#if isSubmitting}
              Please wait; computing your predictions could take several minutes <Spinner
                class="text-center"
                color="blue"
              />
            {:else}
              Submit
            {/if}
          </button>

          {#if message}
            <Alert color={success ? "green" : "red"}>
              <p class="message">{message}</p>
            </Alert>
          {/if}
        </form>
      </div>
    </AccordionItem>
    {#if results}
      <AccordionItem open={results}>
        {#snippet header()}Results{/snippet}
        <div style="display: flex; margin-bottom: 100px;">
          <div style="height:400px;width:400px">
            <Plot grid frame>
              <Line
                data={Object.values(ms2performance)}
                x="recall"
                y="precision"
              />
              <Dot
                data={[ms2performance[`${threshold[0].toFixed(3)}`]]}
                x="recall"
                y="precision"
              />
            </Plot>
          </div>
          <div style="margin-left: 50px; width: 50%;">
            <h3 style="text-align: center; font-weight: bold;">Performance</h3>
            <div style="width: 200px; display: block; margin: 10px auto;">
              <div style="text-align: center;">
                Threshold: {threshold[0].toFixed(3)}
              </div>
              <Slider bind:value={threshold} min="0" max="0.99" step="0.005" />
            </div>

            <div style="display:block; margin: 0 auto; width: 200px;">
              <div>
                <b>Recall:</b>
                {performanceAtThreshold["recall"].toFixed(3)}
              </div>
              <div>
                <b>Precision:</b>
                {performanceAtThreshold["precision"].toFixed(3)}
              </div>
              <div>
                <b>F1-score:</b>
                {(
                  (2 *
                    (performanceAtThreshold["recall"] *
                      performanceAtThreshold["precision"])) /
                  (performanceAtThreshold["recall"] +
                    performanceAtThreshold["precision"])
                ).toFixed(3)}
              </div>
            </div>
          </div>
        </div>
        <div id="table-wrapper">
          <table id="results-table">
            <thead>
              <tr>
                <th>Protein</th>
                <th>Position</th>
                <th>MethylSight 2.0 Score</th>
                <th>Methylated<br />(at 0.75 Precision)</th>
                <th>Methylated<br />(at 0.5 Recall)</th>
                <th>Methylated<br />(at custom threshold)</th>
              </tr>
            </thead>
            <tbody>
              {#each results as r}
                <tr>
                  <td>{r.protein}</td>
                  <td>{r.position}</td>
                  <td>{r.score.toFixed(3)}</td>
                  <td
                    ><span
                      class={r.score > CONSERVATIVE_THRESHOLD ? "green" : "red"}
                      >{r.score > CONSERVATIVE_THRESHOLD ? "✓" : "✗"}</span
                    ></td
                  >
                  <td
                    ><span
                      class={r.score > PERMISSIVE_THRESHOLD ? "green" : "red"}
                      >{r.score > PERMISSIVE_THRESHOLD ? "✓" : "✗"}</span
                    ></td
                  >
                  <td
                    ><span class={r.score > threshold[0] ? "green" : "red"}
                      >{r.score > threshold[0] ? "✓" : "✗"}</span
                    ></td
                  >
                </tr>
              {/each}
            </tbody>
          </table>
        </div>
        <button
          on:click={downloadAsCSV}
          type="submit"
          disabled={isSubmitting}
          class="text-center"
          color="blue"
          style="display: block; margin: 10px auto;"
        >
          Download as CSV
        </button>
      </AccordionItem>
    {/if}
  </Accordion>
  <p class="font-bold text-center">
    You can read a pre-print describing MethylSight 2.0 <a href="">here</a>.
  </p>

  <p>
    If you use MethylSight 2.0 in your work, please cite us (<a
      on:click={() => downloadPlainTextCitation("bib")}>bib</a
    >
    or <a on:click={() => downloadPlainTextCitation("txt")}>plain text</a>).
  </p>

  <p>
    Free predictions of methylation sites within your protein of interest are
    made possible by HuggingFace, which hosts the MethylSight 2.0 model.
  </p>
  <a href="https://huggingface.co/spaces">
    <img id="huggingfacelogo" src="/huggingface.svg" />
  </a>
</main>

<style>
  .container {
    max-width: 900px;
    margin: auto;
    padding: 2rem;
    font-family: system-ui, sans-serif;
  }

  .hero {
    text-align: center;
    margin-bottom: 2rem;
  }

  .hero h1 {
    font-size: 2.5rem;
    color: #2c3e50;
  }

  .hero p {
    font-size: 1.2rem;
    color: #555;
  }

  .upload-form {
    display: flex;
    flex-direction: column;
    gap: 1rem;
  }

  textarea {
    width: 100%;
    font-family: monospace;
    padding: 1rem;
    border: 1px solid #ccc;
    border-radius: 4px;
  }

  input[type="file"] {
    border: none;
  }

  button {
    padding: 0.75rem;
    background-color: #3498db;
    color: white;
    border: none;
    border-radius: 4px;
    cursor: pointer;
    font-weight: bold;
  }

  button[disabled] {
    background-color: #95a5a6;
    cursor: not-allowed;
  }

  .message {
    margin-top: 1rem;
    font-weight: bold;
    color: #2d3436;
  }

  #huggingfacelogo {
    width: 300px;
    display: block;
    margin: 0 auto;
  }

  h1 {
    margin-bottom: 0px;
  }

  h2 {
    font-size: 1.2rem;
  }

  p {
    margin: 1rem;
  }

  .red {
    color: red;
  }

  .green {
    color: green;
  }

  #table-wrapper {
    height: 300px;
    overflow-y: scroll;
  }

  #result-table {
    width: 100%;
  }

  #results-table th {
    width: 20%;
  }

  #results-table td {
    text-align: center;
  }

  #results-table tr {
    width: 100%;
  }

  #results-table thead {
    border-top: 1px solid black;
    border-bottom: 1px solid black;
  }

  #results-table tbody {
    border-bottom: 1px solid black;
  }
</style>
