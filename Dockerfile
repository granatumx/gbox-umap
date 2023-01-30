FROM granatumx/gbox-py-sdk:1.0.0

RUN pip install umap-learn
RUN apt-get install -y libarpack2-dev
RUN pip install scanpy

COPY . .

RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz package.yaml yamls/*.yaml
