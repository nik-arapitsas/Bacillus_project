
for dir in /home/nik_arapitsas/Desktop/test/SRL*/; do
    if [ -d "$dir" ]; then
        shortname=$(basename "$dir")
        mkdir /home/nik_arapitsas/Desktop/test/SRL*/"$shortname"
        echo "Directory /home/nik_arapitsas/Desktop/test/SRL*/$shortname is created"
    fi
done


   