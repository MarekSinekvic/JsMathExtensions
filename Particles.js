class Particle {
	constructor(Position, Velocity, Mass = 1) {
		this.position = Position;
		this.velocity = Velocity;
		this.acceleration = [0,0];
		this.pastPosition = Position;

		this.mass = Mass;

		this.onIntergrateFuncs = [];
	}
	integrate(dt = 1) {
		this.velocity[0] += this.acceleration[0]*dt;
		this.velocity[1] += this.acceleration[1]*dt;

		this.position[0] += this.velocity[0]*dt*1;
		this.position[1] += this.velocity[1]*dt*1;

		// let pastPos = [this.position[0],this.position[1]];
		// this.position[0] = this.position[0]*2 - this.pastPosition[0] + this.acceleration[0]*dt*dt;
		// this.position[1] = this.position[1]*2 - this.pastPosition[1] + this.acceleration[1]*dt*dt;
		// this.pastPosition = pastPos;
		// this.velocity = [this.position[0]-this.pastPosition[0],this.position[1]-this.pastPosition[1]];

		this.acceleration = [0,0];

		for (let func of this.onIntergrateFuncs) {
			func(this);
		}
	}
	addForce(force) {
		this.acceleration[0] += force[0] / this.mass;
		this.acceleration[1] += force[1] / this.mass;
	}
	addIntegrationFunctor(func = () => {}) {
		this.onIntergrateFuncs.push(func);
	}
}
class RigidBody extends Particle {
	constructor(Position, Velocity = [0,0], Rotation = 0, RotationSpeed = 0, Mass = 1, InertiaMoment = 1) {
		super(Position, Velocity, Mass);

		this.position = Position;
		this.velocity = Velocity;
		this.mass = Mass;

		this.rotation = Rotation;
		this.rotSpeed = RotationSpeed;
		this.inertiaMoment = InertiaMoment;

		this.restitution = 0;
		this.friction = 0.5;

		this.rotAcceleration = 0;
		// console.log(this._Energy());
	}
	_Energy() {
		return (this.inertiaMoment*((this.rotSpeed+this.rotAcceleration)**2)+math.magnitude([this.velocity[0]+this.acceleration[0],this.velocity[1]+this.acceleration[1]])**2)/2
	}
	integrate(dt = 1) {
		super.integrate(dt);

		this.rotSpeed += this.rotAcceleration * dt;
		
		this.rotation += this.rotSpeed*dt;// / (Math.PI*2)
		this.rotation %= Math.PI*2;

		this.rotAcceleration = 0;
	}
	getTangentVelocityAt(point) {
		// let M = math.magnitude(point)*0+1;
		return [-point[1] * this.rotSpeed, point[0] * this.rotSpeed];
	}
	getFullVelocityAt(point) {
		let tangentVelocity = this.getTangentVelocityAt(point);
		return [this.velocity[0]+tangentVelocity[0], this.velocity[1]+tangentVelocity[1]];
	}
	addForce(point, force) {
		super.addForce(force);
		let tangent = [-point[1],point[0]];

		// let M = math.magnitude(tangent);
		// if (M == 0) return;
		// console.log(force, tangent, math.dot(force, tangent));
		this.rotAcceleration += math.dot(force, tangent)/this.inertiaMoment;
	}
	addCollisionForce(localPoint, normal, colliderBody, restitution = 0, friction = 0) {
		let normalJ = -this.calcJParameter([this.position[0]+localPoint[0],this.position[1]+localPoint[1]],normal,colliderBody,restitution);
		let tangentJ = -this.calcJParameter([this.position[0]+localPoint[0],this.position[1]+localPoint[1]],[-normal[1],normal[0]],colliderBody,friction);

		this.addForce(localPoint, [normal[0]*normalJ,normal[1]*normalJ]);
		this.addForce(localPoint, [-normal[1]*tangentJ,normal[0]*tangentJ]);
	}
	addStaticCollisionForce(localPoint, normal, restitution = 1, friction = 0, colliderVelocity = [0,0]) {
		this.addCollisionForce(localPoint, normal, {position:[this.position[0]+localPoint[0],this.position[1]+localPoint[1]],mass:10e16,inertiaMoment:10e16,getFullVelocityAt:(r)=>{return colliderVelocity}},restitution,friction);
	}
	calcJParameter(point,normal,colliderBody,restitution = 1) {
		let relativeR1 = [point[0]-this.position[0],point[1]-this.position[1]];
		let relativeR2 = [point[0]-colliderBody.position[0],point[1]-colliderBody.position[1]];

		let relativeV1 = this.getFullVelocityAt(relativeR1);
		let relativeV2 = colliderBody.getFullVelocityAt(relativeR2);
		
		let relativePointVelocity = [relativeV2[0]-relativeV1[0],relativeV2[1]-relativeV1[1]];
		let v1 = -(1+restitution) * math.dot(relativePointVelocity, normal);
		let v2 = (1/this.mass + 1/colliderBody.mass + (relativeR1[0]*normal[1]-relativeR1[1]*normal[0])**2/this.inertiaMoment + (relativeR2[0]*normal[1]-relativeR2[1]*normal[0])**2/colliderBody.inertiaMoment); // r1*n2-n1*r2
		// if (v2 == 0) return 0;
		return v1/v2;
	}
}
class SortedVolumes {
	constructor(XCount = 1, YCount = 1, Bounds = [800,600]) {
		this.scale = [XCount, YCount];
		this.bounds = Bounds;

		this.boxes = [];
		for (let i = 0; i < YCount; i++) {
			this.boxes[i] = [];
			for (let j = 0; j < XCount; j++) {
				this.boxes[i][j] = [];
			}
		}
	}
	IsInBounds(x,y) {
		if (x < 0 || y < 0 || x >= this.scale[0] || y >= this.scale[1]) return false;
		return true;
	}
	GetInVolumePosition(position) {
		if (position[0] < 0 || position[0] > this.bounds[0]) return [0,0];
		if (position[1] < 0 || position[1] > this.bounds[1]) return [0,0];
		return [
			~~(position[0]/this.bounds[0]*this.scale[0]),
			~~(position[1]/this.bounds[1]*this.scale[1])
		];
	}
	AddBodyToVolume(volumePos, body) {
		let ind = this.boxes[volumePos[1]][volumePos[0]].push(body);
		return ind-1;
	}
	RemoveBodyFromVolume(volumePos,name = '') {
		let ind = this.boxes[volumePos[0][1]][volumePos[0][0]].splice(volumePos[1],1);
		for (let i = volumePos[1]; i < this.boxes[volumePos[0][1]][volumePos[0][0]].length; i++) {
			this.boxes[volumePos[0][1]][volumePos[0][0]][i].volumePos[name][1]--;
		}
		return ind;
	}
	TransferInVolume(from, to) {
		if (from[0] == to[0] && from[1] == to[1]) return;
		if (!this.IsInBounds(from[0],from[1])) return;
		if (!this.IsInBounds(to[0],to[1])) return;

		let target = this.boxes[from[1]][from[0]][from[2]];
		if (typeof (target) == 'undefined') console.log(this.boxes[from[1]][from[0]].length, from, to);

		let newIndex = this.boxes[to[1]][to[0]].push(target);
		this.boxes[from[1]][from[0]].splice(from[2],1);


		return newIndex;
	}
	Clear() {
		for (let i = 0; i < this.boxes.length; i++) {
			for (let j = 0; j < this.boxes[i].length; j++) {
				this.boxes[i][j] = [];
			}
		}
	}
	map(func) {
		for (let i = 0; i < this.boxes.length; i++) {
			for (let j = 0; j < this.boxes[i].length; j++) {
				func(this.boxes[i][j],i,j);
			}
		}
	}
	
	AddParticleToVolume(body, name = '') {
		if (typeof (body.volumePos) == 'undefined') body.volumePos = {};
		body.volumePos[name] = [this.GetInVolumePosition(body.position),0];
		body.volumePos[name][1] = this.AddBodyToVolume(body.volumePos[name][0], body);
	}
	OnParticleIntegrate(body,name='') {
		let newVolumePos = this.GetInVolumePosition(body.position);
		// console.log(this.boxes[body.volumePos[name][0][1]][body.volumePos[name][0][0]]);
		if (body.volumePos[name][0][0] != newVolumePos[0] || body.volumePos[name][0][1] != newVolumePos[1]) {
			
			let oldPosition = [body.volumePos[name][0][0],body.volumePos[name][0][1],body.volumePos[name][1]];
			body.volumePos[name][1] = this.TransferInVolume([...body.volumePos[name][0],body.volumePos[name][1]], [...newVolumePos]) - 1;
			body.volumePos[name][0] = [newVolumePos[0],newVolumePos[1]];
			
			let oldVolumeBox = this.boxes[oldPosition[1]][oldPosition[0]];
			for (let i = oldPosition[2]; i < oldVolumeBox.length; i++) {
				// console.log(oldVolumeBox[i].volumePos);
				oldVolumeBox[i].volumePos[name][1]--;
			}
		}
	}
	GetPointNeighbors(point,Radius = 1) {
		let particles = [];
		let cpos = this.GetInVolumePosition(point);
		for (let y = -Radius; y <= Radius; y++) {
			for (let x = -Radius; x <= Radius; x++) {
				if (!this.IsInBounds(cpos[0]+x,cpos[1]+y)) continue;
				
				let box = this.boxes[cpos[1]+y][cpos[0]+x];
				for (let i = 0; i < box.length; i++) {
					particles.push(box[i]);
				}

			}
		}
		return particles;
	}
	GetParticleNeighbors(body, Radius = 1) {
		let particles = [];
		let cpos = this.GetInVolumePosition(body.position);
		for (let y = -Radius; y <= Radius; y++) {
			for (let x = -Radius; x <= Radius; x++) {
				if (!this.IsInBounds(cpos[0]+x,cpos[1]+y)) continue;
				
				let box = this.boxes[cpos[1]+y][cpos[0]+x];
				for (let i = 0; i < box.length; i++) {
					particles.push(box[i]);
				}

			}
		}
		return particles;
	}

	_DebugDraw(colCoff = 10, fill = false) {
		// for (let y = 0; y < this.scale[1]; y++) {
		// 	let p = y * this.bounds[1]/this.scale[1];
		// 	d.line(0, p, this.bounds[0], p, "red");
		// }
		// for (let x = 0; x < this.scale[0]; x++) {
		// 	let p = x * this.bounds[0]/this.scale[0];
		// 	d.line(p, 0, p, this.bounds[1], "red");
		// }
		let s = [this.bounds[0]/this.scale[0],this.bounds[1]/this.scale[1]];
		for (let y = 0; y < this.scale[1]; y++) {
			for (let x = 0; x < this.scale[0]; x++) {
				let v = math.lerp([0,0,255],[255,0,0], this.boxes[y][x].length/colCoff);
				let clr = inRgb(v);
				let px = x * s[0];
				let py = y * s[1];
				if (!fill)
					d.rect(px+1,py+1,s[0]-2,s[1]-2,"black",clr,1,false);
				else
					d.rect(px,py,s[0]-1,s[1]-1,clr,clr,1,true);
			}
		}
		let mPos = this.GetInVolumePosition([input.mouse.x,input.mouse.y]);
		d.txt(this.boxes[mPos[1]][mPos[0]].length.toString(), input.mouse.x,input.mouse.y,"16px Arial","white")
	}
}
class SPHParticle extends Particle {
	constructor(StartPosition,StartVelocity,Radius, Mass = 1) {
		super(StartPosition,StartVelocity, Mass);
		this.radius = Radius;

		this.density = 0;
		this.pressure = 0;
		this.divergence = 0;
		this.heat = 0;

		this.PressureIntensity = 1000;
		this.DefaultDensity = 0.01;
		this.HeatTransferIntensity = 0.001;
		this.Viscosity = 0.1;
		this.Convergence = 1.9;

	}
	Integrate() {
		super.integrate(DT);
	}
	smoothKernel(x) {
		if (x > this.radius) return 0;
		let v1 = (this.radius-x);
		return 4 * (v1*v1*v1)/(this.radius**4);
	}
	smoothKernelDer(x) {
		if (x > this.radius) return 0;
		return 12/(this.radius*this.radius*this.radius*this.radius) * (this.radius-x)*(this.radius-x);
		// return 1;
	}
	smoothKernelDer2(x) {
		if (x > this.radius) return 0;
		return 24/(this.radius*this.radius*this.radius*this.radius) * (this.radius-x);
	}
	addHeat(heat) {
		let oldHeat = this.heat;
		this.heat += heat;
	}
	static EvalIntegrate(particles, volume, DT = 1) {
		let neighbors = [];
		for (let i = 0; i < particles.length; i++) {
			neighbors[i] = volume.GetParticleNeighbors(particles[i],1);//1+Math.floor(particles[i].radius/(volumes.bounds[0]/volumes.scale[0]*1.5))
			particles[i].CalcDensity(neighbors[i]);
		}
		for (let i = 0; i < particles.length; i++) {
			particles[i].CalcForces(neighbors[i]);
		}
		for (let i = 0; i < particles.length; i++) {
			particles[i].CalcDivergence(neighbors[i])
		}
		let pressures = [];
		for (let i = 0; i < particles.length; i++) {
			pressures.push(particles[i].CalcPressure(neighbors[i]));
		}
		for (let i = 0; i < pressures.length; i++) {
			particles[i].pressure = pressures[i];
		}
		for (let i = 0; i < particles.length; i++) {
			particles[i].CalcAdvection(neighbors[i])
		}
	}
	CalcDensity(neighbors) {
		this.density = 0;
		for (let i = 0; i < neighbors.length; i++) {
			if (neighbors[i] == this) continue;
			
			let delta = [neighbors[i].position[0]-this.position[0], neighbors[i].position[1]-this.position[1]];
			let dst = (delta[0]**2+delta[1]**2);
			if (dst < this.radius**2 && dst > 0.01) {
				dst = Math.sqrt(dst);
				this.density += this.smoothKernel(dst, (this.radius+neighbors[i].radius)/2) * neighbors[i].mass;
			}
		}
		this.density *= (1+1*this.heat);
	}
	CalcDivergence(neighbors, DT = 1) {
		this.divergence = 0;
		for (let i = 0; i < neighbors.length; i++) {
			if (neighbors[i] == this || neighbors[i].density < 0.0001) continue;
			
			let delta = [neighbors[i].position[0]-this.position[0], neighbors[i].position[1]-this.position[1]];
			let dst = (delta[0]**2+delta[1]**2);
			if (dst < this.radius**2 && dst > 0.01) {
				dst = Math.sqrt(dst);

				this.divergence += (neighbors[i].mass/neighbors[i].density)*((neighbors[i].velocity[0]-this.velocity[0])*delta[0]/dst+(neighbors[i].velocity[1]-this.velocity[1])*delta[1]/dst)*this.smoothKernel(dst);
			}
		}
	}
	CalcPressure(neighbors) {
		let pressure = 0;
		pressure = (this.density - this.DefaultDensity) * this.PressureIntensity;
		if (pressure < -0) pressure = -0;
		for (let i = 0; i < neighbors.length; i++) {
			if (neighbors[i] == this) continue;

			let delta = [neighbors[i].position[0]-this.position[0], neighbors[i].position[1]-this.position[1]];
			let dst = (delta[0]**2+delta[1]**2);
			if (dst < this.radius**2 && dst > 0.01) {
				dst = Math.sqrt(dst);

				pressure += (neighbors[i].pressure) * this.smoothKernel(dst); //(neighbors[i].mass/neighbors[i].density)*
				let div = this.divergence;
				if (div < -3) div = 3; if (div > 3) div = 3;
				pressure -= this.Convergence * (div); //(neighbors[i].velocity[0]-this.velocity[0] + neighbors[i].velocity[1]-this.velocity[1])
			}
		}
		
		return pressure;
	}
	CalcForces(neighbors) {
		let vd = Math.sqrt(this.velocity[0]**2+this.velocity[1]**2);
		if (vd > 40) this.velocity = [this.velocity[0]/vd*40,this.velocity[1]/vd*40];
		for (let i = 0; i < neighbors.length; i++) {
			if (neighbors[i] == this) continue;
			if (neighbors[i].density < 0.0001 || neighbors[i].density > 10) continue;
			let delta = [neighbors[i].position[0]-this.position[0], neighbors[i].position[1]-this.position[1]];
			let dst = (delta[0]**2+delta[1]**2);
			if (dst < this.radius*this.radius && dst > 0.01) {
				dst = Math.sqrt(dst);

				// let proj = math.dot([neighbors[i].velocity[0]-this.velocity[0], neighbors[i].velocity[1]-this.velocity[1]], [delta[0]/dst,delta[1]/dst]);
				// //if (proj < 0) proj = 0;
				// let Fproj = proj * neighbors[i].mass / neighbors[i].density * this.Viscosity * this.smoothKernelDer2(dst,(this.radius+neighbors[i].radius)/2); 
				
				// this.addForce([delta[0]/dst * Fproj, delta[1]/dst * Fproj]);

				let Fvisc = this.Viscosity * neighbors[i].mass / neighbors[i].density * this.smoothKernelDer2(dst, (this.radius+neighbors[i].radius)/2);
				this.addForce([(neighbors[i].velocity[0]-this.velocity[0])*Fvisc, (neighbors[i].velocity[1]-this.velocity[1])*Fvisc]);

				if (neighbors[i].heat != this.heat)
					this.addHeat(neighbors[i].mass / neighbors[i].density * (neighbors[i].heat-this.heat)/2 * this.smoothKernelDer2(dst, (points[i].radius+neighbors[i].radius)/2) * this.HeatTransferIntensity * DT);
			}
		}
	}
	CalcAdvection(neighbors) {
		for (let i = 0; i < neighbors.length; i++) {
			if (neighbors[i] == this) continue;
			if (neighbors[i].density < 0.0001 || neighbors[i].density > 10) continue;
			let delta = [neighbors[i].position[0]-this.position[0], neighbors[i].position[1]-this.position[1]];
			let dst = (delta[0]**2+delta[1]**2);
			if (dst < this.radius*this.radius && dst > 0.01) {
				dst = Math.sqrt(dst);
				
				let N = (this.pressure+neighbors[i].pressure)/2/ (neighbors[i].density);
				let Fpress = -this.smoothKernelDer(dst,(this.radius+neighbors[i].radius)/2) * neighbors[i].mass * N;
				this.addForce([delta[0]/dst * Fpress,delta[1]/dst * Fpress]);
			}
		}
	}
	
	Bounding(bounds) {
		let normal = [0,0];
		if (this.position[0] > bounds[0]-this.radius)
			normal[0] = -1;
		if (this.position[0] < this.radius)
			normal[0] = 1;
		if (this.position[1] > bounds[1]-this.radius)
			normal[1] = -1;
		if (this.position[1] < this.radius)
			normal[1] = 1;

		if (normal[0] != 0 || normal[1] != 0) {
			if (normal[0] != 0) {
				let p = this.radius+(-normal[0]+1)/2*(bounds[0]-this.radius*2);
				this.position[0] = p;
				this.velocity[0] = -this.velocity[0] * 0.5;
			}
			if (normal[1] != 0) {
				let p = this.radius+(-normal[1]+1)/2*(bounds[1]-this.radius*2);
				this.position[1] = p;
				this.velocity[1] = -this.velocity[1] * 0.5;
			}
			// searchVolume.OnParticleIntegrate(this,'search');
		}
	}
}
class JointParticle extends Particle {
	#spring(target,defaultDistance) {
		let delta = [target.position[0]-this.position[0], target.position[1]-this.position[1]];
		let distance = Math.sqrt(delta[0]**2 + delta[1]**2);
		if (distance == 0) return;
		let direction = [delta[0]/distance,delta[1]/distance];

		let dst = (distance - defaultDistance) / defaultDistance;
		let f = (dst) * this.spring;
		// delta[0]*(1-defaultDistance/distance)

		this.addForce([direction[0]*f, direction[1]*f]);
	}
	#dump(target) {
		let vel = [target.velocity[0]-this.velocity[0],target.velocity[1]-this.velocity[1]];
		let delta = [target.position[0]-this.position[0], target.position[1]-this.position[1]];
		let dst = math.magnitude(delta);
		let dir = [delta[0]/dst,delta[1]/dst];

		let proj = math.dot(vel,dir);

		this.addForce([dir[0]*proj*this.dump,dir[1]*proj*this.dump]);
		return;
		// let delta = [(target.position[0]+target.velocity[0])-(this.position[0]+this.velocity[0]), (target.position[1]+target.velocity[1])-(this.position[1]+this.velocity[1])];
		// if (delta[0] == 0 && delta[1] == 0) return;
		// let M = Math.sqrt(delta[0]**2+delta[1]**2);

		// let D = (M-defaultDistance);
		// // dumpVec[0] += delta[0]/M * D * this.joints[i].dump;
		// // dumpVec[1] += delta[1]/M * D * this.joints[i].dump;
		// this.acceleration[0] += delta[0]/M * D * mass/(1+this.joints.length);
		// this.acceleration[1] += delta[1]/M * D * mass/(1+this.joints.length);
	}
	#wallsIntegrate(dt = 1) {
		const Restitution = 1;
		const Friction = 0.8;

		let normal = [0,0];
		if (this.position[1] <0) {
			this.position[1] = 0;
			normal = [0,1];
		}
		if (this.position[1] > canv.height) {
			this.position[1] = canv.height;
			normal = [0,-1];
		}
		if (this.position[0] < 0) {
			this.position[0] = 0;
			normal = [1,0];
		}
		if (this.position[0] > canv.width) {
			this.position[0] = canv.width;
			normal = [-1,0];
		}
		let tangent = [-normal[1],normal[0]];

		this.acceleration[0] += normal[0] * -math.dot(this.velocity,normal) * Restitution;
		this.acceleration[1] += normal[1] * -math.dot(this.velocity,normal) * Restitution;
		
		this.acceleration[0] += tangent[0] * -math.dot(this.velocity,tangent) * Friction;
		this.acceleration[1] += tangent[1] * -math.dot(this.velocity,tangent) * Friction;
	}
	constructor(Position,Joints = [], Spring = 1, Dump = 0.5) {
		super(Position,[0,0])
		
		this.joints = Joints;
		this.defaultDsts = [];

		this.spring = Spring;
		this.dump = Dump;

		this.isFreezed = false;
	}
	integrate(dt=1) {
		if (this.isFreezed) return;
		super.integrate(dt);
		this.#wallsIntegrate(dt);
	}
	CalcConstraints(dt=1) {
		for (let i = 0; i < this.joints.length; i++) {
			const J = this.joints[i];
			this.#spring(J, this.defaultDsts[i], this.spring, dt);
		}

		for (let i = 0; i < this.joints.length; i++) {
			const target = this.joints[i];
			this.#dump(target, this.dump, dt);
		}
		
	}
	AddJoint(other) {
		let ind = this.joints.push(other);
		this.defaultDsts.push(math.magnitude([other.position[0]-this.position[0], other.position[1]-this.position[1]]));
		return this.joints[ind-1];
	}
	static ExpandMap(map, size = 1) {
		let newMap = [];
		for (let y = 0, ny = 0; y < map.length; y+=1/size, ny++) {
			newMap[ny] = [];
			for (let x = 0, nx = 0; x < map[Math.floor(y)].length; x+=1/size, nx++) {
				newMap[ny][nx] = map[Math.floor(y)][Math.floor(x)];
			}
		}
		return newMap;
	}
	static GenerateSolidObject(map, StartPosition = [0,0], size = 30, jointsRadius = 1, Spring = 0.001, Dump = 0.01) {
		let particles = [];
		let mapLinks = [];
		for (let y = 0; y < map.length; y++) {
			mapLinks[y] = [];
			for (let x = 0; x < map[y].length; x++) {
				if (map[y][x] != 0) {
					let p = new JointParticle([StartPosition[0]+x*size,StartPosition[1]+y*size], [], Spring, Dump);
					particles.push(p);
					mapLinks[y][x] = p;

					if (map[y][x] == 2) {
						p.isFreezed = true;
					}
				} else {
					mapLinks[y][x] = null;
				}
			}
		}
		for (let y = 0; y < map.length; y++) {
			for (let x = 0; x < map[y].length; x++) {
				
				let addJointCount = 0;
				if (map[y][x] == 3) addJointCount = 3; else addJointCount = 0;
				for (let ly = -jointsRadius-addJointCount; ly <= jointsRadius+addJointCount; ly++) {
					for (let lx = -jointsRadius-addJointCount; lx <= jointsRadius+addJointCount; lx++) {
						if (lx == 0 && ly == 0) continue;
						if (y+ly < 0 || y+ly > map.length-1) continue;
						if (x+lx < 0 || x+lx > map[y+ly].length-1) continue;

						// if (ly == -jointsRadius && (lx == -jointsRadius || lx == jointsRadius)) continue;
						// if (ly == jointsRadius && (lx == -jointsRadius || lx == jointsRadius)) continue;

						if (mapLinks[y+ly][x+lx] === null || mapLinks[y][x] === null) continue;


						let j = mapLinks[y][x].AddJoint(mapLinks[y+ly][x+lx]);
						if (map[y+ly][x+lx] == 3) {
							mapLinks[y][x].dump = 6;
							mapLinks[y][x].spring = 0.5;
						}
					}
				}

			}
		}
		
		// let centere = [StartPosition.x+(map[0].length-1)*size/2, StartPosition.y+(map.length-1)*size/2];
		// for (let i = 0; i < newBody.newParticles.length; i++) {
		// 	let delta = [newBody.newParticles[i].position.x - centere[0], newBody.newParticles[i].position.y - centere[1]];
		// 	let M = math.magnitude(delta);
		// 	// newBody.newParticles[i].position.x += 1/(1+Math.abs(delta[1]**1)/200) * 10 * Math.sign(delta[0]);
		// 	newBody.newParticles[i].position.x += (delta[0]/M)  * 0;
		// 	newBody.newParticles[i].position.y += (delta[1]/M)  * 0;
		// }
		
		return particles;
	}
}