# jmfilter
Join Map Locus File filter.

# Installation

```
  git clone https://github.com/theboocock/jmfilter
  cd jmfilter
  python setup.py --user
```
  
# Example - 


``` 
   cd test/ 
   jmfilter -m .7 test_data.txt

```

## Log file should look exactly like

```
2016-05-04 12:51:06,151 INFO Started LOC filter for linkage map creation
2016-05-04 12:51:06,151 INFO Loaded linkage map
2016-05-04 12:51:06,151 INFO Number of loci = 48
2016-05-04 12:51:06,152 INFO Remove data with missing parent information
2016-05-04 12:51:06,152 INFO NLOC before update: 48
2016-05-04 12:51:06,152 INFO NLOC after update: 6
2016-05-04 12:51:06,152 INFO Remove sites that have a bad call rate
2016-05-04 12:51:06,152 INFO NLOC before update: 6
2016-05-04 12:51:06,152 INFO NLOC after update: 6
2016-05-04 12:51:06,153 INFO Look for impossible genotypes
2016-05-04 12:51:06,157 INFO NLOC before update: 6
2016-05-04 12:51:06,157 INFO NLOC after update: 6

```

