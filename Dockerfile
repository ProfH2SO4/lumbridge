FROM python:3.10

ARG APP_NAME=lumbridge

LABEL "image.maintainer"="forgac.matej@gmail.com"
LABEL "sw.maintainer"="forgac.matej@gmail.com"
LABEL "sw.maintainer.website"="https://github.com/ProfH2SO4"

WORKDIR /usr/src/app

COPY . /usr/src/app

RUN pip install --no-cache-dir -r requirements.txt

RUN apt-get update && apt-get install -y \
    wget \
    perl \
    gcc \
    g++ \
    make \
    unzip \
    zip \
    bedtools

# Download and install HOMER
RUN wget http://homer.ucsd.edu/homer/configureHomer.pl \
    && perl configureHomer.pl -install


# Add HOMER bin directory to PATH
ENV PATH="/usr/src/app/bin/:${PATH}"

EXPOSE 80

ENV NAME World

CMD ["python", "./run.py"]
