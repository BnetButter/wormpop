# Use an official Python runtime as a parent image
FROM python:3.10-slim

WORKDIR /app

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir numpy pandas matplotlib scipy

CMD [ "python3", "analysis.py" ]