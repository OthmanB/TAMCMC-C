## What is this

This is Data that were shared to me by Thomas Kallinger on ~ Oct 2023. They are not the one from his paper of 2014, as he lost them.
However, these correspond to a much larger sample of "raw" data following the exact smae kind of analysis as what he carried out in 2014.
This allows to have a reference basis to complement his analysis and to better define how the power spectrum fitting should be carried out.
The model I devised and now implemented in TAMCMC version > 1.84.5 is derived from observation of these data.

## What do you have here

There is basically two files:

- read_kallinger_data.py: Provides means of reading the csv file that Thomas provided to me and to also plot a1, a2 terms as a function of numax for his analysed stars.
- show_kallinger2014.py : This allows you to plot the typical model from Kallinger+2014. It reveals the importance of having 3 harvey instead of two, considering his model.
