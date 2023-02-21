FROM granatumx/gbox-py-sdk:1.0.0

RUN apt-get update
RUN apt-get upgrade -y
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt update
RUN apt-get upgrade -y
RUN apt install -y python3.10
RUN apt install -y python3.10-distutils
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.10 2
RUN curl -sS https://bootstrap.pypa.io/get-pip.py | python3.10

RUN pip install -U pandas
RUN pip install -U numpy
RUN pip install -U numba
RUN pip install -U numexpr
RUN pip install -U pynndescent
RUN pip install -U tensorflow
RUN pip install -U scipy
RUN pip install -U umap-learn
RUN pip install -U pip
RUN pip install -U matplotlib

COPY . .

RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz package.yaml yamls/*.yaml
