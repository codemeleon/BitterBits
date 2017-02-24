# What is it??

This containt the solution of problems I came accross and decided to create a copy for future reference. Most of the information collected from internet (links included).

## Generate list using pandas groupby
```
df = pd.DataFrame( {'a':['A','A','B','B','B','C'], 'b':[1,2,5,5,4,6]})
df.groupby('a')['b'].apply(list)
```

Source: [Stackoverflow](http://stackoverflow.com/questions/22219004/grouping-rows-in-list-in-pandas-groupby)
