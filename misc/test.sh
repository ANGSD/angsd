./realSFS print list1.saf.idx |grep -E "^18\t" >r18
./realSFS print list1.saf.idx -r 18 >en
./realSFS print list1.saf.idx -r 18 >to
cmp en r18
cmp en to

head -n1 r18|md5
./realSFS print list1.saf.idx -r 18:13999898 |md5
./realSFS print2 list1.saf.idx -r 18:13999898 |md5


head -n1 r18|md5
./realSFS print list1.saf.idx -r 18:13999898-13999898 |md5
./realSFS print2 list1.saf.idx -r 18:13999898-13999898 |md5

head -n2 r18|md5
./realSFS print list1.saf.idx -r 18:13999898-13999899 |md5
./realSFS print2 list1.saf.idx -r 18:13999898-13999899 |md5


head -n3 r18|md5
./realSFS print list1.saf.idx -r 18:13999898-13999900 |md5
./realSFS print2 list1.saf.idx -r 18:13999898-13999900 |md5


head -n3 r18|tail -n2|md5
./realSFS print list1.saf.idx -r 18:13999899-13999900 |md5
./realSFS print2 list1.saf.idx -r 18:13999899-13999900 |md5


tail -n3 r18|head -n1|md5
./realSFS print list1.saf.idx -r 18:14100085 |md5
./realSFS print2 list1.saf.idx -r 18:14100085 |md5


tail -n3 r18|head -n1|md5
./realSFS print list1.saf.idx -r 18:14100085-14100085 |md5
./realSFS print2 list1.saf.idx -r 18:14100085-14100085 |md5


tail -n3 r18|head -n2|md5
./realSFS print list1.saf.idx -r 18:14100085-14100086 |md5
./realSFS print2 list1.saf.idx -r 18:14100085-14100086 |md5


tail -n3 r18|head -n3|md5
./realSFS print list1.saf.idx -r 18:14100085-14100087 |md5
./realSFS print2 list1.saf.idx -r 18:14100085-14100087 |md5

tail -n3 r18|head -n3|md5
./realSFS print list1.saf.idx -r 18:14100085-14100088 |md5
./realSFS print2 list1.saf.idx -r 18:14100085-14100088 |md5
