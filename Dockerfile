# Use OpenMC as base image
FROM debian/openmc:latest

# Set working directory
WORKDIR /root

# Copy requirements.txt into the container
COPY requirements.txt .

# Install pip dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Keep container running with bash
CMD ["/bin/bash"]
