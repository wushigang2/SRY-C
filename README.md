# SRY-C
SRY-C is C language version of SRY.  
SRY is my friend wangxiaobo's kmer tool. (https://github.com/caaswxb/SRY)

# Installation
```sh
git clone https://github.com/wushigang2/SRY-C.git
cd SRY-C
make
```

# run SRY-C
```sh
SRY-C
```
```sh
commands:
 zhangya         1 <= k <=  32.
 zy             33 <= k <=  64.
 syn            65 <= k <=  96.
 songyanni      97 <= k <= 128.
```

# run SRY-C zhangya
```sh
SRY-C zhangya
```
```sh
usage: SRY-C zhangya [options]
 -k <int>    kmer size, 1 <= k <= 32, [21]
 -r <string> input kmer file, [NULL]
 -i <string> input read file, [NULL]
 -o <string> output file, [NULL]
```

# run SRY-C songyanni
```sh
SRY-C songyanni
```
```sh
usage: SRY-C songyanni [options]
 -k <int>    kmer size, 97 <= k <= 128, [100]
 -r <string> input kmer file, [NULL]
 -i <string> input read file, [NULL]
 -o <string> output file, [NULL]
```

# Contact
Shigang Wu wushigang@caas.cn  
Xiaobo Wang wangxiaobo@caas.cn
