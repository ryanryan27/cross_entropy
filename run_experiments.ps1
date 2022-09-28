param (
$d = 'd',
$n = 10,
$m = 5,
$a = 0.5,
$s = 0,
$i = 1,
$r = 5,
$h = 0.1
)

$index_offset = 1

$path = ".\graphs\"

$csvout = ".\results\ce_$(get-date -f "yyy-MM-dd_hh-mm-ss").csv"
Add-content $csvout "d=$($d),n=$($n),m=$($m),a=$($a),r=$($r),R=$($h),s=$($s),i=$($i)"
Get-ChildItem -path "$($path)*.txt" | Where-Object {[int]($_.BaseName.split("_")[-3]) -lt 200} | ForEach-Object {$_.BaseName;$op = (&".\cross_entropy.exe" -f "$($path)$($_.Name)" $index_offset -d $d -n $n -m $m -a $a -r $r -R $h -s $s -o -1 -i $i); Add-Content $csvout $op.split("\")[-1]}

