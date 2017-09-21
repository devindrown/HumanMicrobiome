# Running mothur in Batch mode

Mothur is a command line program. You can run it in interactive mode by double-clicking on the mothur.exe executable.

Alternatively, you can run it from the command line.

**To open a command prompt in windows select the Start menu and type cmd**

**Navigate to where the mothur executable lives**

```
cd Documents/mothur
```

It is perfectly acceptable to enter the commands for your analysis from within mothur. We call this the interactive mode. If you are doing a lot these types of analysis without thinking too much there is an alternative.

## Batch mode ##

In the `mothur` folder there is a file called `uafdemodata.batch`

You can open this file in `notepad` and you'll see all of the commands you ran, but instead of listing out the file names it uses the current option throughout. 

The beauty of the batch mode is that you can run mothur from your command line without much typing. 

**To run mothur in batch mode give it the name of the file containing the commands you wish to have executed**
```
mothur.exe uafdemodata.batch
```
Sit back and wait and let it rip. All of the output is going in a folder called `testrunbatch`

The other wonderful thing about this approach is that you can use this very file changing the name of the input sequences.

Once this is finished, it will produce a new file called `uafdemobatch.biom`

**Upload this `.biom` file to phinch.org and explore**
