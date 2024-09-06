from python:latest

WORKDIR /app
COPY . .
RUN pip3 install -r requirements.txt
