pop1True <- scan("pops1.frq")
pop2True <- scan("pops2.frq")

sfs1 <- scan("POPS1.sfs.ml")
sfs2 <- scan("POPS2.sfs.ml")

combTrue <- scan("pops.frq")
combSfs <- scan("POPS.sfs.ml")

pdf("firstRes.pdf")
barplot(rbind(pop1True,sfs1),main="pop1",col=1:2)
legend("topright",c("True","est"),col=1:2,lwd=2)
barplot(rbind(pop2True,sfs2),main="pop2",col=1:2)
legend("topright",c("True","est"),col=1:2,lwd=2)

barplot(rbind(combTrue,combSfs[25]),main="combined (removed cat-25 from dirty.sfs)",col=1:2)
legend("topright",c("True","est"),col=1:2,lwd=2)
dev.off()
