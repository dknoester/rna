source("~/research/src/rna/math/common.R")
setwd("/Users/dk/research/src/rna/var")
figpath = "/Users/dk/research/src/rna/figures/"

D = load.files("multiobjective_fitness.dat.gz", "ta0", "003-nsga2")
summary(D)

g = ggplot(data=D, aes(x=update, y=max_fitness-1)) + stat_summary(aes(color=factor(objective),fill=factor(objective)),fun.data="mean_cl_boot", geom="smooth") + labs(x="Update", y="Fitness (w)") + quick_theme

showfig(g)








D = load.files("rna_multiobjective_detail.dat", "ta0", "003-nsga2")
summary(D)
g = ggplot(data=D, aes(x=update,y=max_fitness,color=objective)) + geom_point() + labs(y="Fitness (w)") + quick_theme
showfig(g)

D = load.files("multiobjective_fitness.dat.gz", "ta0", "003-nsga2")
summary(D)
g = ggplot(data=D, aes(x=update, y=max_fitness-1)) + stat_summary(aes(color=factor(objective),fill=factor(objective)),fun.data="mean_cl_boot", geom="smooth") + labs(x="Update", y="Fitness (w)") + quick_theme
showfig(g)

D[which.max(D$max_fitness),]

plot(m0$max_fitness,m1$max_fitness)
plot(m1$fitness,m2$fitness)
plot(m2$fitness,m3$fitness)
plot(m3$fitness,m4$fitness)

D = load.files("multiobjective_fitness.dat.gz", "ta0", "003-nsga2")
m0 = subset(D,objective==0)
m1 = subset(D,objective==1)
m2 = subset(D,objective==2)
m3 = subset(D,objective==3)
m4 = subset(D,objective==4)

m0[which.max(m0$max_fitness),]
m1[which.max(m1$max_fitness),]
m2[which.max(m2$max_fitness),]
m3[which.max(m3$max_fitness),]
m4[which.max(m4$max_fitness),]

# # > m0[which.max(m0$max_fitness),]
      # update mean_generation objective min_fitness mean_fitness max_fitness                                      filename      expr
# 39991   1895         156.534         0      0.3671       0.6287      0.6584 003-nsga2/ta0_7/multiobjective_fitness.dat.gz 003-nsga2
      # treatment trial
# 39991       ta0     7
# > m1[which.max(m1$max_fitness),]
      # update mean_generation objective min_fitness mean_fitness max_fitness                                       filename      expr
# 21947   1895         954.682         1       0.481       0.9619      0.9661 003-nsga2/ta0_22/multiobjective_fitness.dat.gz 003-nsga2
      # treatment trial
# 21947       ta0    22
# > m2[which.max(m2$max_fitness),]
    # update mean_generation objective min_fitness mean_fitness max_fitness                                      filename      expr treatment
# 978    975         720.392         2      0.6654       0.8544      0.8565 003-nsga2/ta0_1/multiobjective_fitness.dat.gz 003-nsga2       ta0
    # trial
# 978     1
# > m3[which.max(m3$max_fitness),]
     # update mean_generation objective min_fitness mean_fitness max_fitness                                      filename      expr
# 1584   1580        1040.026         3       0.196       0.9668      0.9702 003-nsga2/ta0_1/multiobjective_fitness.dat.gz 003-nsga2
     # treatment trial
# 1584       ta0     1
# > m4[which.max(m4$max_fitness),]
      # update mean_generation objective min_fitness mean_fitness max_fitness                                       filename      expr
# 25925   1860         600.568         4      0.6495       0.8428      0.8463 003-nsga2/ta0_26/multiobjective_fitness.dat.gz 003-nsga2
      # treatment trial
# 25925       ta0    26






g = ggplot(data=D, aes(x=update, y=max_fitness, color=objective)) + geom_point() + quick_theme


stat_summary(aes(color=objective,fill=objective),fun.data="mean_cl_boot", geom="smooth") + labs(x="Update", y="Fitness (w)") + quick_theme
showfig(g)






# Fitness data for the different experiments.
#
#
D = load.files("fitness.dat.gz", "ta0", "004-rnaca")
D = load.files("fitness.dat.gz", "ta0", "005-expca")

g = ggplot(data=D, aes(x=update, y=mean_fitness)) + stat_summary(aes(color=expr,fill=expr),fun.data="mean_cl_boot", geom="smooth") + labs(x="Update", y="Fitness (w)") + quick_theme
showfig(g)

dominant(D)

R = load.files("rna.dat.gz", "ta0", "004-rnaca")
summary(R)
g = ggplot(data=R, aes(x=max_fitness, y=max_test)) + geom_point() + ylim(0.0,1.0) + xlim(0,1.0)
showfig(g)





# Fitness of base experiments.
F_base_fitness = expr_fitness(x004)
# F_base_fitness = F_base_fitness + scale_fill_discrete(breaks=c("002-1d-fsm","003-2d-fsm","004-3d-fsm"),labels=c("1D","2D","3D")) + scale_color_discrete(breaks=c("002-1d-fsm","003-2d-fsm","004-3d-fsm"),labels=c("1D","2D","3D"))
showfig(F_base_fitness)
savefig(F_base_fitness, "p-fitness")


# 1000x dominant fitness
#
#
D = load.files("ca_dom_1000x")
D = subset(D, expr=="002-1d-fsm" | expr=="003-2d-fsm" | expr=="004-3d-fsm")
levels(D$expr) = c("1D", "2D", "3D", "1D-rl", "2D-rl", "3D-rl")
# D$delta_w = D$w1 - D$w0
d = ggplot(D, aes(x=expr,y=w1,fill=expr)) + geom_boxplot() + geom_jitter() + ylim(0,1) + quick_theme + ylab("Fraction correct (out of 1000)") + theme(legend.position="none") + xlab("")
showfig(d)
savefig(d, "p-1000x-fitness")

# D1 = subset(D,expr=="1D")
# D1[which.max(D1$w1),]
# individual     w0    w1                              filename expr treatment trial delta_w
# 6     494817 0.9918 0.865 .//002-1d-fsm/ta0_14/ca_dom_1000x.dat   1D       ta0    14 -0.1268

# D1 = subset(D,expr=="2D")
# D1[which.max(D1$w1),]
# individual     w0    w1                              filename expr treatment trial delta_w
# 52     498839 0.9955 0.884 .//003-2d-fsm/ta0_29/ca_dom_1000x.dat   2D       ta0    29 -0.1115

# D1 = subset(D,expr=="3D")
# D1[which.max(D1$w1),]
# individual     w0    w1                              filename expr treatment trial delta_w
# 63     333564 0.9918 0.881 .//004-3d-fsm/ta0_11/ca_dom_1000x.dat   3D       ta0    11 -0.1108


# Rule table density data
#
#
R = load.files("ca_dom_rule_density")
levels(R$expr) = c("1D")
R = merge(R, subset(S,expr=="1D" & scale==1), by="trial")
# cor(R$w, R$rho, method="spearman")
# == 0.106

P = load.files("ca_dom_sampled_rule_density")
levels(P$expr) = c("2D", "3D")
P = rbind(data.frame(w=P$w, rho=P$rho, expr=P$expr), data.frame(w=R$w, rho=R$rho, expr=R$expr.x))
r = ggplot(P, aes(x=w,y=rho,color=expr)) + geom_point() + quick_theme  + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position="right") + ylab("Rule density") + xlab("Fitness (w)")
showfig(r)
savefig(r, "p-rule-density")

# KO hidden data
K = load.files("ca_dom_ko_hidden.dat")
summary(K)
K$delta_w = K$w1 - K$w0
levels(K$expr) = c("1D", "2D", "3D")
k = ggplot(K, aes(x=expr,y=delta_w,fill=expr)) + geom_boxplot() + geom_jitter() + ylim(-1,0) + quick_theme + ylab("Change in classification accuracy") + theme(legend.position="none") + xlab("")
showfig(k)
savefig(k, "p-ko-hidden")

# KO input data
K = load.files("ca_dom_ko_input.dat")
summary(K)
K$sop = (K$w0 - K$w1) / (K$w0+K$w1)
levels(K$expr) = c("1D", "2D", "3D")
k = ggplot(K, aes(x=expr,y=sop,fill=expr)) + geom_boxplot() + geom_jitter() + ylim(0,.5) + quick_theme + ylab(expression(paste("Self-organization ",(S[op])))) + theme(legend.position="none") + xlab("")
showfig(k)
savefig(k, "p-ko-input")

# Scaling data
#
#
S = load.files("ca_dom_scale.dat", "ta0")
S = subset(S, expr=="002-1d-fsm" | expr=="003-2d-fsm" | expr=="004-3d-fsm")
levels(S$expr) = c("1D", "2D", "3D", "1D-rl", "2D-rl", "3D-rl")
# s = ggplot(S, aes(x=factor(scale),y=w,color=expr)) + geom_point() + labs(x="Scaling factor",y="Fitness") + quick_theme + theme(legend.position="right") + geom_jitter()
# showfig(s)
# savefig(s, "p-scaling-fitness", width=12)


S1 = subset(S,expr=="1D")
S1$cells = 35 * S1$scale
S2 = subset(S,expr=="2D")
S2$cells = (7 * S2$scale) * (7 * S2$scale)
S3 = subset(S,expr=="3D")
S3$cells = (3 * S3$scale) * (3 * S3$scale) * (5 * S3$scale)

s = ggplot(rbind(S1,S2,S3), aes(x=cells,y=w,color=expr)) + geom_point() + geom_jitter() + scale_x_log10() + quick_theme + theme(legend.position="right") + labs(x="Number of cells", y="Fraction correct (out of 1000)")
showfig(s)
savefig(s, "p-scaling-fitness")


# S1 = subset(S, expr=="1D" & scale==9)
# S1[which.max(S1$w),]
# scale    w                              filename expr treatment trial
# 54     9 0.86 .//002-1d-fsm/ta0_14/ca_dom_scale.dat   1D       ta0    14

# S1 = subset(S, expr=="2D" & scale==9)
# S1[which.max(S1$w),]
# scale   w                              filename expr treatment trial
# 441     9 0.9 .//003-2d-fsm/ta0_26/ca_dom_scale.dat   2D       ta0    26

# S1 = subset(S, expr=="3D" & scale==9)
# S1[which.max(S1$w),]
# scale    w                              filename expr treatment trial
# 747     9 0.93 .//004-3d-fsm/ta0_30/ca_dom_scale.dat   3D       ta0    30

