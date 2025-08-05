
# Running MethylSight

## Requirements

- A Unix-style platform (ideally, Linux)
- [Docker](https://docs.docker.com/desktop/setup/install)

## With Docker

1. Install Docker on your system.

2. Build our Docker image (it is large, at ~15 GB).

```sh
git clone https://github.com/GreenCUBIC/MethylSight2.git
cd MethylSight2
docker build -t methylsight2 .
```

3. Create an empty (blank) output file where the results should be stored.

```sh
touch /path/to/my/results/file.csv
```

4. Run MethylSight 2.0

If a GPU is available:

``` sh
input=<ABSOLUTE_PATH_TO_INPUT_FILE>;output=<ABSOLUTE_PATH_TO_OUTPUT_FILE>;docker run -v "$input:/input.fasta:ro" --mount type=bind,source="$output",target="/output.csv" --gpus 1 methylsight2 "/env/bin/python3 model.py -i /input.fasta -w weights.ckpt -o /output.csv"
```

If no GPUs are available:
``` sh
input=<ABSOLUTE_PATH_TO_INPUT_FILE>;output=<ABSOLUTE_PATH_TO_OUTPUT_FILE>;docker run -v "$input:/input.fasta:ro" --mount type=bind,source="$output",target="/output.csv" methylsight2 "/env/bin/python3 model.py -i /input.fasta -w weights.ckpt -o /output.csv"
```

The model weights can be downloaded [here](https://methylsight2.s3.us-east-2.amazonaws.com/weights.ckpt).

## Huggingface

``` sh
curl -X POST -H "Content-Type: application/json" -d '{"sequence": "SEQUENCE"}' https://fcharih-methylsight2.hf.space/predict
```
