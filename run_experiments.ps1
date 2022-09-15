param (
$d = 't',
$n = 10,
$m = 5,
$a = 0.5,
$s = 0,
$i = 1
)

$index_offset = 1

$path = ".\graphs\"

$csvout = ".\results\ce_$(get-date -f "yyy-MM-dd_hh-mm-ss").csv"

Get-ChildItem -path "$($path)*.txt" | ForEach-Object {$op = (&".\cross_entropy.exe" -f "$($path)$($_.Name)" $index_offset -d $d -n $n -m $m -a $a -s $s -o -1 -i $i); Add-Content $csvout $op.split("\")[-1]}

