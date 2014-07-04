# Solution based on: http://unix.stackexchange.com/a/25144
# But two improvements:
# 1. added header row 
# 2. print columns in user-specified order as per http://stackoverflow.com/a/10666771/3806692
BEGIN {
    n=split(cols,out,",")
    for (i=1; i<=n; i++)
        printf "%s%s", out[i], OFS
    print ""
}
NR==1 {
    for (i=1; i<=NF; i++)
        ix[$i] = i
}
NR>1 {
    for (i=1; i<=n; i++)
        printf "%s%s", $ix[out[i]], OFS
    print ""
}
