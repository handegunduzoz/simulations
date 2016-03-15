

y
 nopr
 prop,,1



 tecp,,-1
 tecp


plot,pers,1
0
0 0 300
0,1,0

loop,,1
dt,,1.d-1
time
loop,,20
tang,,1
next
plot,hide
plot,cont,5
 tecp
save
next


loop,,16
dt,,1.d0
time
loop,,20
tang,,1
next
 tecp,,5
plot,wipe
plot,hide
plot,cont,5
save
next


dt,,1.d-3
loop,,50
time
loop,,20
tang,,1
next
 tecp,,5
save
plot,wipe
plot,hide
plot,cont,5
next


dt,,1.d-4
loop,,50
time
loop,,20
tang,,1
next
 tecp,,5
save
plot,wipe
plot,hide
plot,cont,5
next

dt,,1.d-5
loop,,50
time
loop,,20
tang,,1
next
 tecp,,5
save
plot,wipe
plot,hide
plot,cont,5
next



