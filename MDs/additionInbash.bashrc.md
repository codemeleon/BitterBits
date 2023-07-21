i# Addded by Anmol. Please do not change it

if [ ${USER} != "root" ]
then
export PATH="/.anmol/anaconda3/bin:$PATH"
alias ariba="/.anmol/anaconda3/envs/py365/bin/ariba"
alias srst2="/.anmol/anaconda3/envs/py27/bin/srst2"

fi


echo "


Welcome ${USER}

Always backup your data and source codes, especially the results and source code you need for publications and presentations.

If you want admin to backup up your codes automatically. Please put all your code in '${USER}' folder in your home directory and execute code from there. All the code will be uploaded to github.com/codemeleon/${USER} public repository. Please also learn a pipeline tools such as Ruffus [http://www.ruffus.org.uk/], CWL [https://www.commonwl.org/], NextFlow [https://www.nextflow.io/] for reproducible results and easy code sharing.  Always create a README.md file for  each project and update it regularly.

If you have any query, Please contact Anmol Kiran (akiran@mlw.mw).

Thank you for your cooperation.



"

