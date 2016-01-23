awk '!/^(│|.─|Exiting)/ {gsub(/sage\: /, ""); print}'  output | awk 'BEGIN{ RS="OK"; FS="\n"}; {split($2,a, " "); split($6, b, " ") ; print a[3] "," b[2] }' | awk '/[0-9]+,[0-9]+/ { print  }'
