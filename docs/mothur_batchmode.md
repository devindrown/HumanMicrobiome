# Running mothur in Batch mode

It is perfectly acceptable to enter the commands for your analysis from within mothur. We call this the interactive mode. If you are doing a lot these types of analysis without thinking too much there is an alternative.

## Batch mode ##

In the `Home` folder there is a file called `mothur_MiSeq_SOP.batch`

You can open this file by clicking on the name and you'll see all of the commands you ran, but instead of listing out the file names it uses the current option throughout. 

The beauty of the batch mode is that you can run mothur from your command line without much typing. 

**To run mothur in batch mode give it the name of the file containing the commands you wish to have executed**
```
mothur mothur_MiSeq_SOP.batch
```
Sit back and wait and let it rip. All of the output is going in a folder called `testrun_batch`

The other wonderful thing about this approach is that you can use this very file changing the name of the input sequences.

Once this is finished, it will produce a new file called `testrun_BATCH.biom`

**Load this `.biom` file to phinch.org and explore**
