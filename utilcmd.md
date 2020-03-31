### Parse .csv files where fields contain commas

```
perl -MText::ParseWords -lne '@a = parse_line(",", 0, $_); print join("\t", map { $_ } @a)' input.csv >output.tab
```

