FROM python:3.10

WORKDIR /usr/src/app

# Copy the current directory contents into the container at /usr/src/app
COPY . /usr/src/app

RUN pip install --no-cache-dir -r requirements.txt

# Update and install necessary packages
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

# Define environment variable
ENV NAME World

CMD ["python", "./run.py"]
