## TOOLS to be added here

You will need to add JAFFA here
git clone https://github.com/Oshlack/JAFFA.git

Now edit the following variables in JAFFA/JAFFA_stages.groovy
```
codeBase = "/absolute/path/to/JAFFA"
refBase = "\$GENOMES/GRCh37/jaffa"
genome = "GRCh37"
annotation = "75"
genomeFasta = "\$GENOMES/GRCh37/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa"
jaffa_output = "jaffa/"
```

Ensure all tools are executable.
```
find . -type f -exec chmod 700 '{}' \;
```
