Performance tests were conducted on computers equipped with a 2.10GHz Intel(R) Xeon(R) Gold 6230R CPU and 32GB of RAM in simulated WAN (with 200Mbps bandwidth and 20ms latency) and LAN (with 2500Mbps bandwidth, 1ms latency) environments.


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

# Run 2PC-based protocols


- "lh" is the length of the histogram, or M in our paper. 
- "nd" is the number of data points, or "N" in our paper. 
- "nt" is the number of threads.
- The bit width is automatically set to a value that is sufficient for the given parameters. It can also be set manually using the "bw" parameter.

``` bash
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

```
