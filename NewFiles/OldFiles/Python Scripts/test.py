h = 2
dt = 2

init_vel = 1
init_acc = 2
b_values = 3


vel1 = init_vel + ((h*dt)*init_acc) + ((h*h*dt/2)*b_values) + ((h*h*h*dt/3)*b_values) + ((h*h*h*h*dt/4)*b_values)+ ((h*h*h*h*h*dt/5)*b_values)+ ((h*h*h*h*h*h*dt/6)*b_values)+ ((h*h*h*h*h*h*h*dt/7)*b_values)+ ((h*h*h*h*h*h*h*h*dt/8)*b_values)+ ((h*h*h*h*h*h*h*h*h*dt/9)*b_values);
	
vel2 = init_vel+ ((h*dt)   *	(init_acc+ ((h/2)    * 	(b_values+ ((2*h/3)  * 	(b_values + ((3*h/4)  * 	(b_values+ ((4*h/5) *	(b_values+ ((5*h/6)*	(b_values+ ((6*h/7)*	(b_values+ ((7*h/8)*	(b_values+ ((8*h/9)*	(b_values))))))))))))))))));	


