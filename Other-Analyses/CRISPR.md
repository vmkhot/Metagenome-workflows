# Analyses of CRISPR Arrays

### minced

### CRISPRCasfinder

### [CRISPRviz](https://github.com/CRISPRlab/CRISPRviz)

There is some ssh-gymnastics involved in getting this to run and opening the results from the server because of a firewall. Use mobaxterm. Open an ssh tunnel with 4444 ports on both side. when this is established, localhost:4444 will open results in your own browser

```bash
# requested resources to go into a node/shell
 salloc --mem=10G -c 5 -N 1 -t 05:00:00
 ```

 Ran crisprviz on node for 1 fasta file but it went through all fasta files in the directory?? It produces a lot of files, each crispr goes into its own repeat and spacer file. This command also froze the shell but all the things were output.

```bash
crisprviz.sh -x -f bin3_crispr_fullcontig.fasta
```

To open, I ran this command so that it would produce the html output

```bash
crisprviz.sh -s bin3_crispr_fullcontig.fasta_spacers.fa
```

Note: this will only work if you have a ssh tunnel between remote server and your local computer.  

Then in your browser: http://localhost:4444 to open the results. The spacers are very messy?? I don't know what to do with these visuals and you can only download a screenshot - very inconvenient.
