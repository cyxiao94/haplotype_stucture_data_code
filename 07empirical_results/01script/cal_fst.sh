#split
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,4-18 > 04fst/anc_F0.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,19-33 > 04fst/23_F10.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,34-48 > 04fst/23_F20.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,49-63 > 04fst/23_F30.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,64-78 > 04fst/23_F40.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,79-93 > 04fst/23_F50.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,94-108 > 04fst/23_F60.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,109-123 > 04fst/23_F70.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,124-138 > 04fst/23_F80.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,139-153 > 04fst/2818_F10.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,154-168 > 04fst/2818_F20.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,169-183 > 04fst/2818_F30.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,184-198 > 04fst/2818_F40.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,199-213 > 04fst/2818_F50.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,214-228 > 04fst/2818_F60.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,229-243 > 04fst/2818_F70.sync
zcat 00data/max15.mix85.SA.hfluc_hcons.F0-F80.flt.sync.gz | cut -f 1-3,244-258 > 04fst/2818_F80.sync

input=$1

#anc
perl ~/Changyi/softwares/popoolation2_1201/fst-sliding.pl --input 04fst/anc_F0.sync --output 04fst/anc_F0.fst --min-count 2 --min-coverage 4 --max-coverage 2% --pool-size 600 --window-size 1000 --step-size 1000
rm 04fst/anc_F0.sync

#23
for i in 10 20 30 40 50 60 70 80
do
  perl ~/Changyi/softwares/popoolation2_1201/fst-sliding.pl --input 04fst/23_F${i}.sync --output 04fst/23_F${i}.fst --min-count 2 --min-coverage 4 --max-coverage 2% --pool-size 600 --window-size 1000 --step-size 1000   
  rm 04fst/23_F${i}.sync
done

#2818
for i in 10 20 30 40 50 60 70 80
do
  perl ~/Changyi/softwares/popoolation2_1201/fst-sliding.pl --input 04fst/2818_F${i}.sync --output 04fst/2818_F${i}.fst --min-count 2 --min-coverage 4 --max-coverage 2% --pool-size 600 --window-size 1000 --step-size 1000   
  rm 04fst/2818_F${i}.sync
done    