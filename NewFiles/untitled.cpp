pos[j] = init_pos(i,j)
	+ ((h*dt) * init_vel(i,j))
	+ ((h*h*dt*dt/2)*init_acc(i,j)) 
	+ ((h*h*h*dt*dt/6)*b_values[i](0,j)) 
	+ ((h*h*h*h*dt*dt/12)*b_values[i](1,j)) 
	+ ((h*h*h*h*h*dt*dt/20)*b_values[i](2,j))
	+ ((h*h*h*h*h*h*dt*dt/30)*b_values[i](3,j))
	+ ((h*h*h*h*h*h*h*dt*dt/42)*b_values[i](4,j))
	+ ((h*h*h*h*h*h*h*h*dt*dt/56)*b_values[i](5,j))
	+ ((h*h*h*h*h*h*h*h*h*dt*dt/72)*b_values[i](6,j))
	+ ((h*h*h*h*h*h*h*h*h*h*dt*dt/90)*b_values[i](7,j));
	
vel[j] = init_vel(i,j)
	+ ((h*dt)   *	(init_acc(i,j)
	+ ((h/2)    * 	(b_values[i](0,j)
	+ ((2*h/3)  * 	(b_values[i](1,j) 
	+ ((6*h/4)  * 	(b_values[i](2,j)
	+ ((24*h/5) *	(b_values[i](3,j)
	+ ((120*h/6)*	(b_values[i](4,j)
	+ ((720*h/7)*	(b_values[i](5,j)
	+ ((5040*h/8)*	(b_values[i](6,j)
	+ ((40320*h/9)*	(b_values[i](7,j)))))))))))))))))));	


