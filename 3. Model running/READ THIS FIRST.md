First, load the workspace saved after 1. 

Then, run plotClarkeGrid_adj.R first. 

Next, run any of the other model-running scripts. LR and DT include personal + population models, XGBoost and RF are run separately. 
After running a script, save the workspace again, and reload the workspace saved after 1. This gives a relatively clean workspace, since especially RF is very memory intensive and likely will crash R. 

Metrics and feature importance plots are automatically saved. I do recommend creating save states at the end of the for loop. 
With XGBoost and RF I also recommend adding save states after each training. This makes it easier to stop the code, or if R crashes, you still have the trainings.

To investigate the best hyperparameters, check the best_model at the end, this should provide the information
