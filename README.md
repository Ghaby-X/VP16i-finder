# VP16i-finder

## A web interface for predictive models of VP16 inhibitors built using flask

This project provides an interface for preictive models that has been developed to predict the inhibition activity of VP16. It takes in a Canonical smile as input and returns an activity class and the confidence score of the predictive models.

The models were built on the dataset provided by PubChem AID 651615.


## How to run the project locally
To run this project, you would have to have python installed and follow the steps below;

1. Clone the project in your preferred directory
   ```
   git clone 
   ```
3. 
4. Navigate to project folder
5. activate the virtual environment
for cmd or powershell
code -> venv\Scripts\activate
6. Install the libraries in requirement.txt file
pip install -r requirements.txt
7. Run the following code to start the project
flask run
8. Go to the link showing on your terminal

## Usage
for single molecule prediction,
1. Go to upload page <- the first section of upload page
2. Enter a valid smile
3. Click on submit

for multiple molecule prediction;
1. Go to upload page <- Second section of upload page
2. Click on choose file and upload a csv containing your canonical smiles. Your csv must conform to the structure provided in example.csv
3. Click on submit

## Visual representation of the project
