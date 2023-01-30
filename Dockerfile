FROM granatumx/gbox-py-sdk:1.0.0

RUN apt-get update
RUN apt-get upgrade -y

RUN pip install -U umap-learn scanpy python-igraph leidenalg scipy
RUN pip install -U numba
RUN pip install -U numpy
RUN pip install -U numexpr

COPY . .

RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz package.yaml yamls/*.yaml
