Performance tests were conducted on computers equipped with a 2.10GHz Intel(R) Xeon(R) Gold 6230R CPU and 32GB of RAM in simulated WAN (with 200Mbps bandwidth and 20ms latency) and LAN (with 2500Mbps bandwidth, 1ms latency) environments.

> HE相关的在myseal项目中的ndss.cpp和U-test.cpp中，加解密性能用官方的sealbench测试。MBFV用lattigo测试。

## LAN

Setting

``` bash
sudo tc qdisc add dev ens18 root netem delay 0.7ms rate 2500mbit
```

Delete

``` bash
sudo tc qdisc del dev ens18 root netem delay 0.7ms rate 2500mbit
```

## WAN

Setting

``` bash
sudo tc qdisc add dev ens18 root netem delay 20ms rate 200mbit
```

Delete

``` bash
sudo tc qdisc del dev ens18 root netem delay 20ms rate 200mbit
```

# 主要命令

``` bash
# "lh" is the length of the histogram, or M in our paper. "nd" is the number of data points, or "N" in our paper. "nt" is the number of threads.
# The bit width is automatically set to a value that is sufficient for the given parameters. It can also be set manually using the "bw" parameter.
./PSAC_OT_quantile_hideK r=1 ip="{host2's ip}" lh=4000 nd=4000 nt=1 #host1
./PSAC_OT_quantile_hideK r=2 ip="{host1's ip}" lh=4000 nd=4000 nt=1 #host2

./PSAC_GC_quantile_hideK r=1 ip="{host2's ip}" lh=4000 nd=4000 nt=1 sgc=64 #host1
./PSAC_GC_quantile_hideK r=2 ip="{host1's ip}" lh=4000 nd=4000 nt=1 sgc=64 #host2

./PSAC_OT_quantile_rmind r=1 ip="{host2's ip}" lh=4000 nd=4000 nt=1 #host1
./PSAC_OT_quantile_rmind r=2 ip="{host1's ip}" lh=4000 nd=4000 nt=1 #host2

./PSAC_GC_sort_rmind r=1 ip="{host2's ip}" lh=4000 nd=4000 sgc=16
./PSAC_GC_sort_rmind r=2 ip="{host1's ip}" lh=4000 nd=4000 sgc=16

./PSAC_OT_outlier r=1 ip="{host2's ip}" lh=4000 nd=10000 nt=4
./PSAC_OT_outlier r=2 ip="{host1's ip}" lh=4000 nd=10000 nt=4


./my-U-OT r=1 ip="{host2's ip}" lh=4000 nd=100 nt=4
./my-U-OT r=2 ip="{host1's ip}" lh=4000 nd=100 nt=4

./ndss O=4000 N=100

# ./PSAC_GC_chi2 r=1 ip="192.168.8.117" nd=4000 M=2 L=50 sgc=16
# ./PSAC_GC_chi2 r=2 ip="192.168.8.180" nd=4000 M=2 L=50 sgc=16


```

