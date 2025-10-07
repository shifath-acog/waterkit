# Use the existing image as a base to inherit all scientific dependencies
FROM waterkit-cli

# Set the working directory inside the container
WORKDIR /home/waterkit/WATERKIT

# Copy and install the new API-specific requirements
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Copy the API code into the container
# Note: We copy the whole directory to match the python module path
COPY . /home/waterkit/WATERKIT/api

# Command to run the API server
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
