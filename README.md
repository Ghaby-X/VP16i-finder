# VP16i-finder

## A web interface for predictive models of VP16 inhibitors built using flask

This project provides an interface for preictive models that has been developed to predict the inhibition activity of VP16. VP16 is a transactivator viral protein in the the herpes simplex various which contributes to the transcription of immediate early genes in the virus lifecycle. This project takes in a Canonical smile as input and returns an activity class and the confidence score of the predictive models.

The models were built on the dataset provided by PubChem AID 651615.


## How to run the project locally
To run this project, you would have to have python installed and follow the steps below;

1. Clone the project in your preferred directory
   
   ```
   git clone https://github.com/Ghaby-X/VP16i-finder.git
   ```
2. Navigae to the project folder

   ```
   cd VP16i-finder
   ```
   
3. activate virtual environment

   ```
   venv\Scripts\activate
   ```

6. Install the libraries in requirement.txt file
   ```
   pip install -r requirements.txt
   ```
   
8. Run the following code to start the project
   
   ```
   flask run
   ```
   
9. Go to the link showing on your terminal

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

https://github.com/Ghaby-X/VP16i-finder/assets/105595126/69500b84-95f0-4c85-986f-bc9db1a32a74


## Project Challenges
The models were built as prediction for inhibitors of VP16, however in the experiment performed, the protein used was a chimeric protein of GAL4-BDB and VP16 TAD. This implies the specificity of the models in the prediction of VP16 would be low.

