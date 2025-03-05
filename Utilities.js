var math = new (class _ {
	constructor() {
		for (let i = 0; i < 10000; i++) {
			this.randoms.push(Math.random());
		}
	}

	normalize(x) {
		if (Array.isArray(x)) {
			let magn = this.magnitude(x);
			if (magn == 0) return [0,0];
			let y = [];
			for (let i = 0; i < x.length; i++) {
				y[i] = x[i] / magn;
			}
			return y;
		}
		if (typeof (x.z) == "number") {
			let d = Math.sqrt(x.x * x.x + x.y * x.y + x.z * x.z);
			if (d == 0) return [0,0];
			return {
				x: x.x / d,
				y: x.y / d,
				z: x.z / d
			};
		}
		let d = Math.sqrt(Math.pow(x.x, 2) + Math.pow(x.y, 2));
		if (d == 0) return [0,0];
		return {
			x: x.x / d,
			y: x.y / d
		};
	}

	sqrMagnitude(x) {
		return Math.pow(x.x, 2) + Math.pow(x.y, 2);
	}
	magnitude(x) {
		if (Array.isArray(x)) {
			let l = 0;
			for (let i = 0; i < x.length; i++) {
				l += x[i] * x[i];
			}
			return Math.sqrt(l);
		} else {
			let d = Math.pow(x.x, 2) + Math.pow(x.y, 2);
			if (x.z != null) d += Math.pow(x.z, 2);
			return Math.sqrt(d);
		}
	}
	dot(a, b) {
		if (Array.isArray(a)) {
			let v = 0;
			for (let i = 0; i < a.length; i++) {
				v += a[i] * b[i];
			}
			return v;
		}
		if (typeof (a.z) == "number") {
			if (typeof (a.w) == "number")
				return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
			return a.x * b.x + a.y * b.y + a.z * b.z;
		}
		return a.x * b.x + a.y * b.y;
	}
	cross(a, b) {
		// let M1 = new Matrix([
		// 					[0,-a.z,a.y],
		// 					[a.z,0,-a.x],
		// 					[-a.y,a.x,0]]);
		// let M2 = new Matrix([[b.x], [b.y], [b.z]]);

		let c = {
			x: -a.z * b.y + a.y * b.z,
			y: a.z * b.x - a.x * b.z,
			z: -a.y * b.x + a.x * b.y
		};

		return c;
	}
	pseudoDot(a, b) {
		if (Array.isArray(a)) {
			return a[0]*b[1]-b[0]*a[1];
		}
		return a.x * b.y - b.x * a.y;
	}
	reflect(v, normal, K = 2) {
		normal = math.normalize(normal);

		let dotProduct = v.x * normal.x + v.y * normal.y;
		dotProduct = Math.abs(dotProduct);

		return { x: v.x - K * dotProduct * normal.x, y: v.y - K * dotProduct * normal.y };
	}
	distanceToLine(target, p1, p2) {
		if (p1.length == 2) {
			let delta = [p2[0] - p1[0], p2[1] - p1[1]];
			let M = Math.sqrt(delta[0] ** 2 + delta[1] ** 2);
			// if (M == 0) return Math.sqrt((target[0] - p1[0]) ** 2 + (target[1] - p1[1]) ** 2 + (target[2] - p1[2]) ** 2);
			let N = [delta[0] / M, delta[1] / M];

			let tInLine = -(N[0] * (p1[0] - target[0]) + N[1] * (p1[1] - target[1]));
			// tInLine /= M;

			if (tInLine < 0) {
				// return math.magnitude([p1[0] - target[0], p1[1] - target[1]]);
				return [Math.sqrt((p1[0] - target[0]) ** 2 + (p1[1] - target[1]) ** 2),tInLine];
			}
			if (tInLine > M) {
				// return math.magnitude([p2[0] - target[0], p2[1] - target[1]]);
				return [Math.sqrt((p2[0] - target[0]) ** 2 + (p2[1] - target[1]) ** 2),tInLine];
			}
			// let tFromLine = math.dot(N, [-target[1] + p1[1], target[0] - p1[0]]);
			// let tFromLine = N[0] * (p1[1] - target[1]) + N[1] * (target[0] - p1[0]);
			let tFromLine = Math.sqrt(((target[0] - p1[0]) ** 2 + (target[1] - p1[1]) ** 2) - tInLine ** 2);
			return [Math.abs(tFromLine),tInLine];
		} else if(p1.length == 3) {
			let delta = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
			let M = Math.sqrt(delta[0] ** 2 + delta[1] ** 2 + delta[2] ** 2);
			// if (M == 0) return Math.sqrt((target[0] - p1[0]) ** 2 + (target[1] - p1[1]) ** 2 + (target[2] - p1[2]) ** 2);
			let N = [delta[0] / M, delta[1] / M, delta[2] / M];

			let tInLine = -(N[0] * (p1[0] - target[0]) + N[1] * (p1[1] - target[1]) + N[2] * (p1[2] - target[2]));
			// tInLine /= M;

			if (tInLine < 0) {
				// return math.magnitude([p1[0] - target[0], p1[1] - target[1]]);
				return [Math.sqrt((p1[0] - target[0]) ** 2 + (p1[1] - target[1]) ** 2 + (p1[2] - target[2]) ** 2), tInLine];
			}
			if (tInLine > M) {
				// return math.magnitude([p2[0] - target[0], p2[1] - target[1]]);
				return [Math.sqrt((p2[0] - target[0]) ** 2 + (p2[1] - target[1]) ** 2 + (p2[2] - target[2]) ** 2), tInLine];
			}
			// let tFromLine = math.dot(N, [-target[1] + p1[1], target[0] - p1[0]]);
			// let tFromLine = N[0] * (p1[1] - target[1]) + N[1] * (target[0] - p1[0]);
			let tFromLine = Math.sqrt(((target[0] - p1[0]) ** 2 + (target[1] - p1[1]) ** 2 + (target[2] - p1[2]) ** 2) - tInLine ** 2);
			return [Math.abs(tFromLine),tInLine];
		}
		return "dimension"
	}
	linesIntersectionData(P1,P2,p1,p2) {
		let mainLineLength = math.magnitude([P2[0]-P1[0],P2[1]-P1[1]]);
		let delta = [p2[0]-p1[0],p2[1]-p1[1]];
		let targetLineLength = math.magnitude(delta);
		let T = math.pseudoDot(delta,[p1[0]-P1[0],p1[1]-P1[1]])/math.pseudoDot(delta,[P2[0]-P1[0],P2[1]-P1[1]]);

		let intersectPoint = [P1[0]+(P2[0]-P1[0])*T,P1[1]+(P2[1]-P1[1])*T];
        let inLineT = math.dot(math.normalize(delta),[intersectPoint[0]-p1[0],intersectPoint[1]-p1[1]]);

		let isIntersect = true;
		if ((T < 0 || T > mainLineLength) || (inLineT < 0 || inLineT > targetLineLength)) isIntersect = false;

		return {T: T, isIntersect: isIntersect, inLineT: inLineT};
	}

	Triangulate(P, test = P.length) {
		let eachPointConnectedLines = [];
		for (let i = 0; i < P.length; i++) {
			eachPointConnectedLines[i] = [];
		}
		let lines = [];
		for (let ind = 0; ind < test; ind++) {
			let targetPoint = P[ind];
			let targetPointIndex = ind;
			let lockedPoints = [];

			// while (true) {
				let min1 = Infinity,min2 = Infinity;
				let min1Index = -1, min2Index = -1;
				for (let i = 0; i < P.length; i++) {
					if (i==targetPointIndex) continue;

					let delta = [P[i][0]-targetPoint[0],P[i][1]-targetPoint[1]];
					let distance = (delta[0]**2) + (delta[1]**2);
					
					let haveIntersection = false;
					for (let j = 0; j < lines.length; j++) {
						let intersection;
						if (targetPointIndex != lines[j][0] && i != lines[j][1])
							intersection = this.linesIntersectionData(targetPoint,P[i],P[lines[j][0]],P[lines[j][1]]);
						else
							continue;
						if (intersection.isIntersect) {
							haveIntersection = true;
							break;
						}
					}

					if (!haveIntersection) {
						if (distance < min1) {
							min2 = min1;
							min2Index = min1Index;

							min1 = distance;
							min1Index = i;
						} else if (distance < min2) {
							min2 = distance;
							min2Index = i;
						}

					}
				}
				if (min1 != Infinity && min2 != Infinity) {
					lines.push([targetPointIndex,min1Index]);
					lines.push([targetPointIndex,min2Index]);
					lines.push([min1Index,min2Index]);

					eachPointConnectedLines[targetPointIndex].push(min1Index);
					eachPointConnectedLines[targetPointIndex].push(min2Index);
					eachPointConnectedLines[min1Index].push(min2Index);

					eachPointConnectedLines[min1Index].push(targetPointIndex);
					eachPointConnectedLines[min2Index].push(targetPointIndex);
					eachPointConnectedLines[min2Index].push(min1Index);
				} else {
					lockedPoints.push(targetPointIndex);
				}
			// }
		}
		return lines;
	}
	sort(arr, f) {
		let transTable = [];
		let newArray = [];
		for (let i = 0; i < arr.length; i++) {
			newArray[i] = arr[i];
			transTable[i] = i;
		}
		let i = 0;

		let state = 1;
		while (i < arr.length - 1) {
			if (f(newArray[i], newArray[i + 1])) {
				let e1 = newArray[i];
				newArray[i] = newArray[i + 1];
				newArray[i + 1] = e1;

				let e2 = transTable[i];
				transTable[i] = transTable[i + 1];
				transTable[i + 1] = e2;

				i = 0;
				continue;
			}
			i++;
		}
		return {sorted: newArray, transTable: transTable};
	}
	Triangulate2(P, steps = P.length) {
		let lines = [];
		for (let i = 0; i < steps; i++) {
			console.log(i +": ");
			console.log("Lines count: " + lines.length);
			let PSortedByDistance = math.sort(P,(a,b)=>{
				let dst1 = math.magnitude([a[0]-P[i][0],a[1]-P[i][1]]);
				let dst2 = math.magnitude([b[0]-P[i][0],b[1]-P[i][1]]);
				// console.log(P);
				return dst1 > dst2;
			});
			PSortedByDistance.sorted.splice(0,1);
			PSortedByDistance.transTable.splice(0,1);


			// for (let j = 1; j < PSortedByDistance.sorted.length; j++) {
			// 	let isIntersectSome = false;
			// 	for (let k = 0; k < lines.length; k++) {
			// 		if ((i == lines[k][0] || i == lines[k][1]) && (PSortedByDistance.transTable[j] == lines[k][0] || PSortedByDistance.transTable[j] == lines[k][1])) {
			// 			console.log("when j = " + j + " ("+PSortedByDistance.transTable[j]+") and k = " + k + " in equal");
			// 			isIntersectSome = true;
			// 			break;
			// 		}
			// 		if (i == lines[k][0] || i == lines[k][1]) continue;
			// 		if (PSortedByDistance.transTable[j] == lines[k][0] || PSortedByDistance.transTable[j] == lines[k][1]) continue;

			// 		let intersection1 = math.linesIntersectionData(P[i],PSortedByDistance.sorted[j],P[lines[k][0]],P[lines[k][1]]);
			// 		let intersection2 = math.linesIntersectionData(PSortedByDistance.sorted[j-1],PSortedByDistance.sorted[j],P[lines[k][0]],P[lines[k][1]]);

			// 		// if (PSortedByDistance.transTable[0] == lines[k][0] || PSortedByDistance.transTable[0] == lines[k][1]) continue;
			// 		// if (PSortedByDistance.transTable[1] == lines[k][0] || PSortedByDistance.transTable[1] == lines[k][1]) continue;
			// 		if ((intersection1.isIntersect && intersection1.T != 0) || (intersection2.isIntersect && intersection2.T != 0)) {
			// 			isIntersectSome = true;
			// 			console.group();
			// 			console.log("when j = " + j + " ("+PSortedByDistance.transTable[j]+") it intersect with " + k);
			// 			console.log(intersection1);
			// 			console.log(intersection2);
			// 			console.groupEnd();
			// 			// console.log(P[i],PSortedByDistance.sorted[j],P[lines[k][0]],P[lines[k][1]]);
			// 			break;
			// 		}
			// 	}
			// 	if (!isIntersectSome) {
			// 		lines.push([ i, PSortedByDistance.transTable[j-1] ]);
			// 		lines.push([ i, PSortedByDistance.transTable[j] ]);
			// 		lines.push([ PSortedByDistance.transTable[j-1], PSortedByDistance.transTable[j] ]);
			// 	}
			// }
			let isIntersectSome = false;
			for (let k = 0; k < lines.length; k++) {
				// if (lines[j][0] == i || lines[j][1] == i)

				if (lines[k][0] == PSortedByDistance.transTable[1] || lines[k][1] == PSortedByDistance.transTable[1]) continue;
				if (lines[k][0] == PSortedByDistance.transTable[0] || lines[k][1] == PSortedByDistance.transTable[0]) continue;

				let intersection1 = math.linesIntersectionData(P[i],PSortedByDistance.sorted[0],P[lines[k][0]],P[lines[k][1]]);
				let intersection2 = math.linesIntersectionData(P[i],PSortedByDistance.sorted[1],P[lines[k][0]],P[lines[k][1]]);
				let intersection3 = math.linesIntersectionData(PSortedByDistance.sorted[0],PSortedByDistance.sorted[1],P[lines[k][0]],P[lines[k][1]]);

				console.group("Data for k = " + k);
				console.log(intersection1);
	 			console.log(intersection2);
				if ((intersection1.isIntersect) || (intersection2.isIntersect) || intersection3.isIntersect) {
					isIntersectSome = true;
			 			// console.log("when j = " + 1 + " ("+PSortedByDistance.transTable[1]+") it intersect with " + k);
			 			// console.log(intersection1);
			 			// console.log(intersection2);
					// console.log(P[i],PSortedByDistance.sorted[j],P[lines[k][0]],P[lines[k][1]]);
 					console.groupEnd();
					break;
				} else {
				}
 				console.groupEnd();
			}
			if (!isIntersectSome) {
				lines.push([ i, PSortedByDistance.transTable[0] ]);
				lines.push([ i, PSortedByDistance.transTable[1] ]);
				lines.push([ PSortedByDistance.transTable[0], PSortedByDistance.transTable[1] ]);
			}
		}
		return lines;
	}
	///-BASIC-///

	floor(x) {
		return Math.floor(x);
	}
	floorWithNegate(x) {
		return Math.sign(x) * Math.floor(Math.abs(x));
	}
	round(x) {
		return Math.round(x);
	}
	ceil(x) {
		return Math.ceil(x);
	}
	pow(x, y) {
		return Math.pow(x, y);
	}
	sqrt(x) {
		return Math.sqrt(x);
	}
	clamp(x, a, b) {
		if (x < a) return a;
		if (x > b) return b;
		return x;
	}

	///-BASIC-///

	distance(p0, p1) {
		return [Math.sqrt(Math.pow(p0[0] - p1[0], 2) + Math.pow(p0[1] - p1[1], 2)), p0[0] - p1[0], p0[1] - p1[1]];
	}
	lerp(a, b, t) {
		if (Array.isArray(a) && Array.isArray(b)) {
			let c = [];
			for (let i = 0; i < a.length; i++) {
				c.push(a[i] + (b[i] - a[i]) * t);
			}
			return c;
		}
		if (typeof a == "object") {
			return {
				x: a.x + (b.x - a.x) * t,
				y: a.y + (b.y - a.y) * t
			};
		}
		return a + (b - a) * t;
	}

	circleCollision(origin, direction, circlePosition, circleRadius = 1) {
		direction = math.normalize(direction);
		let deltaPosition = {
			x: circlePosition.x - origin.x,
			y: circlePosition.y - origin.y
		};
		let b = 2 * math.dot(direction, deltaPosition);
		let c = math.dot(deltaPosition, deltaPosition) - (circleRadius * circleRadius);
		let D = b * b - 4 * c;
		if (D > 0) {
			let t1 = (b - Math.sqrt(D)) / 2;
			let t2 = (b + Math.sqrt(D)) / 2;
			if (t1 > 0 && t2 > 0) {
				let x1 = {
					x: origin.x + direction.x * t1,
					y: origin.y + direction.y * t1
				};
				let x2 = {
					x: origin.x + direction.x * t2,
					y: origin.y + direction.y * t2
				};
				return [true, x1, x2, t1, t2];
			} else if (t1 < 0 && t2 < 0) {
				return [false];
			} else {
				let x = {
					x: origin.x + direction.x * t2,
					y: origin.y + direction.y * t2
				};
				return [true, x, 0, t2, 0];
			}
		} else if (D == 0) {
			let t = b / 2;
			if (t > 0) {
				let x = {
					x: origin.x + direction.x * t,
					y: origin.y + direction.y * t
				};
				return [true, x, 0, t, 0];
			} else
				return [false];
		} else {
			return [false];
		}
	}
	dstToLine(lx1, ly1, lx2, ly2, x, y) {
		let delta1 = {
			x: lx2 - lx1,
			y: ly2 - ly1
		};
		let normal = math.normalize({ x: -delta1.y, y: delta1.x });
		let delta2 = {
			x: x - lx1,
			y: y - ly1
		};
		return { byNormal: math.dot(delta2, normal), byLine: math.dot(delta2, math.normalize(delta1)), byLineNormalized: (math.dot(delta2, math.normalize(delta1)) / math.magnitude(delta1)) };
	}
	getEqDeriative(eq, x, d = 0.000001) {
		return (eq(x + d) - eq(x)) / d;
	}
	getXzEqDeriative(eq, x, z, d = 0.000001) {
		return {
			x: (eq(x + d, z) - eq(x, z)) / d,
			z: (eq(x, z + d) - eq(x, z)) / d
		};
	}
	frac(a) {
		return a - Math.floor(a);
	}
	// random(x = 0,y = 0,z = 0,w= 0) {
	// 	return this.frac(Math.sin(math.dot({x:x,y:y,z:z,w:w},{x:102.9898,y:708.233,z:153.8465,w:9845.8465}))*43758.5453123);
	// }
	randoms = [];
	getRandomFromArray(i) {
		if (i > 30000) {
			i = i % this.randoms.length;
		}
		if (i > this.randoms.length - 1) {
			let l = this.randoms.length
			for (let j = 0; j < i - l + 1; j++) {
				this.randoms.push(Math.random());
			}
			return this.randoms[i];
		} else {
			return this.randoms[i];
		}
	}
	hash(i) {
		return (i * 651959 + 19698) % 1682;
	}
	random(point) { // point is in N dimensions
		let dot = 0;
		for (let i = 0; i < point.length; i++) {
			dot += point[i] * this.hash(i);
		}
		return math.frac(Math.sin(dot) * 1.);
	}
	RandomVector(DimsCount = 2, seed = -1) {
		let v = [];
		let leng = 0;
		for (let i = 0; i < DimsCount; i++) {
			let c = Math.random();
			v.push(c);
			leng += c * c;
		}
		leng = Math.sqrt(leng);
		for (let i = 0; i < DimsCount; i++) {
			v[i] = v[i] / leng;
		}
		return v;
	}
	SimplexSmooth(x) {
		return x * x * x * (x * (x * 6.0 - 15.0) + 10.0);
	}
	SimplexNoise0(x, y, z = 0, w = 0) { // dims = []
		let pointInGrid = {
			x: Math.floor(x),
			y: Math.floor(y),
			z: Math.floor(z),
			w: Math.floor(w)
		};
		let localPoint = {
			x: x - pointInGrid.x,
			y: y - pointInGrid.y,
			z: z - pointInGrid.z,
			w: w - pointInGrid.w
		};

		let cubeVertexCount = Math.pow(2, 4);
		let cubePoints = [];
		for (let i = 0; i < cubeVertexCount; i++) {
			cubePoints.push({ // 1,2,4,8 = 2 ** n; n - dimension
				x: Math.floor(i / 1) % 2,
				y: Math.floor(i / 2) % 2,
				z: Math.floor(i / 4) % 2,
				w: Math.floor(i / 8) % 2
			});
		}

		let randomVectors = [];
		for (let i = 0; i < cubeVertexCount; i++) {
			let rx = -1 + math.random([(cubePoints[i].x + pointInGrid.x) * 1000, (cubePoints[i].y + pointInGrid.y) * 100, (cubePoints[i].z + pointInGrid.z) * 10, (cubePoints[i].w + pointInGrid.w)]) * 2;
			let ry = -1 + math.random([(cubePoints[i].x + pointInGrid.x) * 100, (cubePoints[i].y + pointInGrid.y) * 1000, (cubePoints[i].z + pointInGrid.z) * 10, (cubePoints[i].w + pointInGrid.w)]) * 2;
			let rz = -1 + math.random([(cubePoints[i].x + pointInGrid.x) * 10, (cubePoints[i].y + pointInGrid.y) * 100, (cubePoints[i].z + pointInGrid.z) * 1000, (cubePoints[i].w + pointInGrid.w)]) * 2;
			let rw = -1 + math.random([(cubePoints[i].x + pointInGrid.x), (cubePoints[i].y + pointInGrid.y) * 10, (cubePoints[i].z + pointInGrid.z) * 100, (cubePoints[i].w + pointInGrid.w) * 1000]) * 2;

			let leng = Math.sqrt(rx * rx + ry * ry + rz * rz + rw * rw);
			randomVectors.push({
				x: rx / leng,
				y: ry / leng,
				z: rz / leng,
				w: rw / leng
			});
		}

		let dots = [];
		for (let i = 0; i < cubeVertexCount; i++) {
			let dir = {
				x: localPoint.x - cubePoints[i].x,
				y: localPoint.y - cubePoints[i].y,
				z: localPoint.z - cubePoints[i].z,
				w: localPoint.w - cubePoints[i].w,
			};
			dots.push(dir.x * randomVectors[i].x + dir.y * randomVectors[i].y + dir.z * randomVectors[i].z + dir.w * randomVectors[i].w);
		}

		localPoint.x = this.SimplexSmooth(localPoint.x);
		localPoint.y = this.SimplexSmooth(localPoint.y);
		localPoint.z = this.SimplexSmooth(localPoint.z);
		localPoint.w = this.SimplexSmooth(localPoint.w);

		let xlerp1 = math.lerp(dots[0], dots[1], localPoint.x);
		let xlerp2 = math.lerp(dots[2], dots[3], localPoint.x);
		let xlerp3 = math.lerp(dots[4], dots[5], localPoint.x);
		let xlerp4 = math.lerp(dots[6], dots[7], localPoint.x);
		let xlerp5 = math.lerp(dots[8], dots[9], localPoint.x);
		let xlerp6 = math.lerp(dots[10], dots[11], localPoint.x);
		let xlerp7 = math.lerp(dots[12], dots[13], localPoint.x);
		let xlerp8 = math.lerp(dots[14], dots[15], localPoint.x);

		let ylerp1 = math.lerp(xlerp1, xlerp2, localPoint.y);
		let ylerp2 = math.lerp(xlerp3, xlerp4, localPoint.y);
		let ylerp3 = math.lerp(xlerp5, xlerp6, localPoint.y);
		let ylerp4 = math.lerp(xlerp7, xlerp8, localPoint.y);

		let zlerp1 = math.lerp(ylerp1, ylerp2, localPoint.z);
		let zlerp2 = math.lerp(ylerp3, ylerp4, localPoint.z);

		let lerp1 = math.lerp(zlerp1, zlerp2, localPoint.w);
		return [lerp1];
	}
	VoronoiNoise0(x,y,z=0,w=0) {
		let pointInGrid = [
			Math.floor(x),
			Math.floor(y),
			Math.floor(z),
			Math.floor(w)
		];
		let localPoint = [
			x - pointInGrid[0],
			y - pointInGrid[1],
			z - pointInGrid[2],
			w - pointInGrid[3]
		];

		let v = 1;
		for (let w = -1; w <= 1; w++) {
			for (let z = -1; z <= 1; z++) {
				for (let y = -1; y <= 1; y++) {
					for (let x = -1; x <= 1; x++) {
						let point = [
							math.random([(x+pointInGrid[0])+1000, (y+pointInGrid[1]), (z+pointInGrid[2]), (w+pointInGrid[3])]),
							math.random([(x+pointInGrid[0]), (y+pointInGrid[1])+1000, (z+pointInGrid[2]), (w+pointInGrid[3])]),
							math.random([(x+pointInGrid[0]), (y+pointInGrid[1]), (z+pointInGrid[2])+1000, (w+pointInGrid[3])]),
							math.random([(x+pointInGrid[0]), (y+pointInGrid[1]), (z+pointInGrid[2]), (w+pointInGrid[3])+1000]),
						];
						let delta = [
							localPoint[0] - (point[0]+x), localPoint[1] - (point[1]+y), localPoint[2] - (point[2]+z), localPoint[3] - (point[3]+w),
						];
						let len = math.magnitude(delta);
						if (v > len) {
							v = len;
						}
					}
				}
			}
		}
		return 1-v;
	}
	SimplexNoise(point) { // dims = []
		let pointInGrid = [];
		for (let i = 0; i < point.length; i++) {
			pointInGrid.push(Math.floor(point[i]));
		}
		let localPoint = [];
		for (let i = 0; i < pointInGrid.length; i++) {
			localPoint.push(point[i] - pointInGrid[i]);
		}

		let cubeVertexCount = Math.pow(2, point.length);
		let cubePoints = [];
		for (let i = 0; i < cubeVertexCount; i++) {
			let cubePoint = [];
			for (let j = 0; j < point.length; j++) {
				cubePoint.push(Math.floor(i / Math.pow(2, j)) % 2);
			}
			cubePoints.push(cubePoint);
		}


		let randomVectors = [];
		for (let i = 0; i < cubeVertexCount; i++) {
			let randomVector = [];
			let leng = 0;
			let randIndex = 0;
			for (let j = 0; j < point.length; j++) {
				let x = (cubePoints[i][j] + pointInGrid[j]);
				randIndex += Math.pow(10, j) * x;
			}
			for (let j = 0; j < point.length; j++) {
				let x = (cubePoints[i][j] + pointInGrid[j]);
				// let v = -1 + this.random([randIndex]) * 2;
				let v = -1 + this.getRandomFromArray(Math.abs(randIndex)) * 2; // randIndex**2
				randomVector.push(v);
				leng += v * v;
			}
			leng = Math.sqrt(leng);
			for (let j = 0; j < randomVector.length; j++) {
				randomVector[j] /= leng;
			}
			randomVectors.push(randomVector);
		}

		let dots = [];
		for (let i = 0; i < cubeVertexCount; i++) {
			let dir = [];
			for (let j = 0; j < point.length; j++) {
				dir.push(localPoint[j] - cubePoints[i][j]);
			}
			let dot = 0;
			for (let j = 0; j < point.length; j++) {
				dot += dir[j] * randomVectors[i][j];
			}
			dots.push(dot);
		}

		for (let i = 0; i < localPoint.length; i++) {
			localPoint[i] = this.SimplexSmooth(localPoint[i]);
		}

		let lerps = [];
		for (let i = 0; i < Math.pow(2, point.length); i += 2) {
			lerps.push(this.lerp(dots[i], dots[i + 1], localPoint[0]));
		}

		for (let i = 0; lerps.length > 1; i++) {
			let newLerp = [];
			for (let j = 0; j < lerps.length; j += 2) {
				newLerp.push(math.lerp(lerps[j], lerps[j + 1], localPoint[i + 1]));
			}
			lerps = newLerp;
		}
		return lerps[0];
	}
	FBM(point, octavesCount = 1) {
		let v = 0;
		let freq = 1, intens = 1;
		for (let t = 0; t < octavesCount; t++) {
			v += this.SimplexNoise(this.mult(point, freq)) * intens;
			freq *= 2;
			intens *= 0.5;
		}
		return v / 1;
	}
	mult(a, b) {
		if (Array.isArray(a) && typeof (b) == "number") {
			for (let i = 0; i < a.length; i++) {
				a[i] *= b;
			}
			return a;
		}
		return a;
	}
	isArrayHaveChar(char, arr) {
		for (let i = 0; i < arr.length; i++) {
			if (arr[i] == char) {
				return true;
			}
		}
		return false;
	}
})();

class Vector {
	values;
	offset;
	length;
	#normalizeArgsType(...args) {
		if (Array.isArray(args[0])) {
			let v2 = [];
			if (typeof (args[1]) != 'undefined')
				v2 = args[1];
			else
				v2 = new Array(args[0].length).fill(0);
			return [args[0], v2];
		} else if (typeof (args[0]) == 'object') {
			let values = [];
			if (args[0].constructor.name == "Vector") {
				for (let i = 0; i < args[0].length; i++) {
					values[i] = args[0][i];
				}
			} else {
				if (args[0].hasOwnProperty('x'))
					values[0] = args[0].x;
				if (args[0].hasOwnProperty('y'))
					values[1] = args[0].y;
				if (args[0].hasOwnProperty('z'))
					values[2] = args[0].z;
				if (args[0].hasOwnProperty('w'))
					values[3] = args[0].w;
			}
			
			let v2 = [];
			if (typeof (args[1]) != 'undefined')
				v2 = args[1];
			else
				v2 = new Array(args[0].length).fill(0);

			return [values, v2];
		} else {
			return [args, new Array(args.length).fill(0)];
		}
	}
	constructor(...args) {
		let needDefineOffset = true;
		this.offset = [];

		let nArgs = this.#normalizeArgsType(...args);
		this.values = nArgs[0];
		this.offset = nArgs[1];
		this.length = this.values.length;

		for (let i = 0; i < this.values.length; i++) {
			this.__defineGetter__(i,()=>{return (this.values[i]+this.offset[i])});
			this.__defineSetter__(i,(v)=>{this.values[i] = v;});
			// this.__defineSetter__();
			// this[i] = this.values[i];
		}
	}
	get() {
		return this.values.map((v,i)=>{return this.values[i]+this.offset[i]});
	}
	set(...args) {
		
	}
	add(v) {
		for (let i = 0; i < v.values.length; i++) {
			this.values[i] += v.values[i];
		}
		return this;
	}
	Draw(origin) {
		
		d.ray(origin[0], origin[1], {x:this[0],y:this[1]}, math.magnitude(this.get()), "red");
	}
}
function replaceInString(str, start, end, target) {
	return str.slice(0,start)+target+str.slice(end+1);
}
function CalcAABB(points) {
	let min = [Infinity,Infinity];
	let max = [-Infinity,-Infinity];
	for (let i = 0; i < points.length; i++) {
		if (points[i][0] < min[0]) min[0] = points[i][0];
		if (points[i][1] < min[1]) min[1] = points[i][1];
		if (points[i][0] > max[0]) max[0] = points[i][0];
		if (points[i][1] > max[1]) max[1] = points[i][1];
	}
	return [min,max];
}
function inRgb(r, g, b, a = 1) {
	if (Array.isArray(r)) {
		let A = 1;
		if (typeof r[3] != "undefined")
			A = r[3];
		return "rgba(" + r[0] + "," + r[1] + "," + r[2] + "," + A + ")";
	}
	return "rgba(" + r + "," + g + "," + b + "," + a + ")";
}
function rgbStringToArray(color) {
	let rgbcolors = color.slice(4, color.length - 1).split(',');
	return [parseFloat(rgbcolors[0]), parseFloat(rgbcolors[1]), parseFloat(rgbcolors[2]), 1]
}
function rgbaStringToArray(color) {
	let rgbcolors = color.slice(5, color.length - 1).split(',');
	return [parseFloat(rgbcolors[0]), parseFloat(rgbcolors[1]), parseFloat(rgbcolors[2]), parseFloat(rgbcolors[3])]
}
function makeShader(gl, src, type) {
    var shader = gl.createShader(type);
    gl.shaderSource(shader, src);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        alert("Error compiling shader: " + gl.getShaderInfoLog(shader));
        return;
    }
    return shader;
}
function initShaders(gl, vs_source, fs_source) {
    // Compile shaders
    var vertexShader = makeShader(gl, vs_source, gl.VERTEX_SHADER);
    var fragmentShader = makeShader(gl, fs_source, gl.FRAGMENT_SHADER);

    // Create program
    var glProgram = gl.createProgram();

    // Attach and link shaders to the program
    gl.attachShader(glProgram, vertexShader);
    gl.attachShader(glProgram, fragmentShader);
    gl.linkProgram(glProgram);
    if (!gl.getProgramParameter(glProgram, gl.LINK_STATUS)) {
        alert("Unable to initialize the shader program");
        return false;
    }

    // Use program
    gl.useProgram(glProgram);
    gl.program = glProgram;

    return true;
}
function initVertexBuffers(gl) {
    // Vertices
    var dim = 2;
    var vertices = new Float32Array([
        -1, 1, 1, 1, 1, -1, // Triangle 1
        -1, 1, 1, -1, -1, -1 // Triangle 2 
    ]);

    // Fragment color

    // Create a buffer object
    var vertexBuffer = gl.createBuffer();
    if (!vertexBuffer) {
        console.log('Failed to create the buffer object');
        return -1;
    }
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.STATIC_DRAW);

    var a_Position = gl.getAttribLocation(gl.program, 'a_Position');
    if (a_Position < 0) {
        console.log('Failed to get the storage location of a_Position');
        return -1;
    }
    gl.vertexAttribPointer(a_Position, dim, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(a_Position);

    gl.uniform2fv(gl.getUniformLocation(gl.program, 'u_Resolution'), [canv.width, canv.height]);
    gl.uniform1i(gl.getUniformLocation(gl.program, 'u_Frame'), frame);

    return vertices.length / dim;
}