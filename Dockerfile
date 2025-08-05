FROM nvidia/cuda:12.9.1-cudnn-runtime-ubuntu24.04 

RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y python3.12 wget
COPY requirements.txt /
RUN apt-get install -y python3-pip
RUN apt install -y python3.12-venv
RUN python3 -m venv env
RUN apt-get install -y pkg-config
RUN apt-get install -y cmake
RUN env/bin/pip install flair
RUN env/bin/pip install -r requirements.txt

RUN mkdir /methylsight2
WORKDIR /methylsight2
ADD https://methylsight2.s3.us-east-2.amazonaws.com/weights.ckpt /methylsight2/methylsight2.ckpt

COPY model.py /methylsight2
COPY server.py /methylsight2

ENTRYPOINT [ "/bin/bash", "-l", "-c" ]