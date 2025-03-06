class TrainableParams3 {
	GetInSortedSpace(coords) {
		let fcoords = [...coords];
		let flatInd = 0;
		for (let i = 0; i < coords.length; i++) {
			if (coords[i]+this.bounds/2 < 0 || coords[i]+this.bounds/2 > this.bounds) continue;
			fcoords[i] = ~~((coords[i]+this.bounds/2)*this.CoverArea);
			
			flatInd += fcoords[i]*((this.bounds*this.CoverArea)**i);
		}

		return ~~(flatInd);
	}
	GetNeighbors(coords) {
		const target = this.GetInSortedSpace(coords);
		let neighs = [...this.coverSortedSpace[target]];
		for (let i = 0; i < Object.keys(this.params[0]).length; i++) {
			// console.log(target,this.bounds*this.CoverArea,i,~~((this.bounds*this.CoverArea)**i));
			// console.log(this.coverSortedSpace[target-~~((this.bounds*this.CoverArea)**i)]);
			const I1 = target+(~~((this.bounds*this.CoverArea)**i));
			const I2 = target-(~~((this.bounds*this.CoverArea)**i));

			if (I1 < 0 || I1 >= this.coverSortedSpace.length) continue;
			if (I2 < 0 || I2 >= this.coverSortedSpace.length) continue;

			neighs = [...neighs,...(this.coverSortedSpace[I1])];
			neighs = [...neighs,...(this.coverSortedSpace[I2])];
		}
		return neighs;
	}
	GetCoverDivergence(point) {
		let localCovers = this.GetNeighbors(point);
		const start = this.CoverSpace(point);
		let div = new Array(point.length).fill(0);
		let sumDst = 0;
		for (const cover of localCovers) {
			let dst = 0;
			for (let i = 0; i < cover[0].length; i++) {
				const delta = cover[0][i]-point[i];
				dst += delta**2;
			}
			sumDst += Math.sqrt(dst);
		}
		// console.log(sumDst);
		// console.log(localCovers);
		
		for (const cover of localCovers) {
			let dst = 0;
			for (let i = 0; i < cover[0].length; i++) {dst += (cover[0][i]-point[i])**2}
			dst = Math.sqrt(dst);
			for (let i = 0; i < cover[0].length; i++) {
				div[i] += (this.CoverSpace(cover[0])-start) * (cover[0][i]-point[i]) / dst; //
			}
		}
		return div;
	}
	constructor(params, scoreFunc, scoreDerFunc, dt = 0.1, coverArea = 0.5) {
		this.DT=dt;
		this.CoverArea = coverArea;
		this.MinDeriative = 10e-10;
		this.MinVelocity = 10e-4;
		this.AgentsCount = 20;
		this.RandomPlaceDispersion = 10;
		this.bounds = 80;

		this.pastParamsDelta = [];
		this.paramsDelta = [];

		this.pastParams = [];
		this.params = [];

		this.state = [];
		for (let i = 0; i < this.AgentsCount; i++) {
			this.pastParams[i] = {};
			this.params[i] = {};
			for (let param in params) {
				let R = (Math.random()-0.5)*2*5;
				this.pastParams[i][param] = params[param]+R;
				this.params[i][param] = params[param]+R;
			}
			this.state[i] = 1;
		}
		this.cover = [];
		this.coverSortedSpace = [];
		for (let i = 0; i < (this.bounds*coverArea)**Object.keys(this.params[0]).length; i++) {
			this.coverSortedSpace.push([]);
		}

		this.onExtrema = false;

		this.oldScore = 0;

		this.scoreFunc = scoreFunc;
		this.scoreDerFunc = scoreDerFunc;
	}
	AddCover(coords,coveredScore,score) {
		if (Math.abs(coveredScore) > 0.05) {
			const nind = this.cover.push([coords,coveredScore,score]) - 1;
			this.coverSortedSpace[this.GetInSortedSpace(coords)].push(this.cover[nind]);
			
			this.cover.sort((a,b)=>{return (b[1])-(a[1])});
		}
	}
	CoverBasis(x) {
		if (Math.abs(x) > 1) return 0;
		// return 1-Math.abs(x);
		return 1/(1+(3.14*x)**2);
	}
	CoverSpace(target) {
		let y = 0;
		// for (let i = 0; i < this.cover.length; i++) {
		// 	let dst = 0;
		// 	for (let j = 0; j < target.length; j++) {
		// 		dst += (target[j] - this.cover[i][0][j])**2;
		// 	}
		// 	dst = Math.sqrt(dst);
		// 	y += this.cover[i][1] * this.CoverBasis(this.CoverArea*dst);
		// }
		let localCovers = this.GetNeighbors(target);
		// console.log(localCovers);
		
		for (let i = 0; i < localCovers.length; i++) {
			let dst = 0;
			for (let j = 0; j < target.length; j++) {
				dst += (target[j] - localCovers[i][0][j])**2;
			}
			dst = Math.sqrt(dst);
			y += localCovers[i][1] * this.CoverBasis(this.CoverArea*dst);
		}
		return y;
	}
	PlaceOnRandomCover(targetIndex) {
		if (this.cover.length == 0) return;

		const gaussDistr = (x,a=1)=>{return Math.exp(-a*Math.PI*(x**2))}
		let targetCover = Math.floor(gaussDistr(Math.random(),this.RandomPlaceDispersion)*(this.cover.length-1));
		if (this.state[targetIndex] == -1) targetCover = this.cover.length-1 - targetCover;
		// console.log(this.cover[0][2]);
		
		
		let ind = 0;
		for (let param in this.params[targetIndex]) {this.params[targetIndex][param] = this.cover[targetCover][0][ind]+(-1+Math.random()*2)*1; ind++}
	}
	IntegrateParams(der, targetIndex) {
		// console.log(der);
		let derMagn = 0;
		for (let param in der) {
			derMagn += der[param];
		}
		
		let index = 0;
		for (let param in this.params[targetIndex]) {
			
			if (Math.abs(derMagn)**0.5 > this.MinDeriative) {
				this.params[targetIndex][param] += der[index] * this.DT * this.state[targetIndex];
			} else {
				// this.PlaceOnRandomCover(targetIndex);
			}
			
			index++;
		}
		
	}
	DeviateParams() {
		for (let i = 0; i < this.AgentsCount; i++) {
			const score = this.scoreFunc(Object.values(this.params[i]));
			const coveredScore = score - this.CoverSpace(Object.values(this.params[i]));
			const der = this.scoreDerFunc(Object.values(this.params[i]));
			this.IntegrateParams(der,i);

			let velMag = 0, isEqual = true;
			this.pastParamsDelta[i] = {}; this.paramsDelta[i] = {};
			for (let param in this.params[i]) {this.pastParamsDelta[i][param] = this.paramsDelta[i][param];}
			for (let param in this.params[i]) {this.paramsDelta[i][param] = this.params[i][param]-this.pastParams[i][param];}
			for (let param in this.params[i]) {velMag += this.paramsDelta[i][param]**2; if (isEqual && this.paramsDelta[i][param] != this.pastParamsDelta[i][param]) isEqual = false;}
			velMag = Math.sqrt(velMag);

			if (isEqual) {
				this.PlaceOnRandomCover(i);
				console.log('equal');
			}
			if (velMag < this.MinVelocity) {
				if (!this.onExtrema) {
					let coords = Object.values(this.params[i]);
					this.AddCover(coords,coveredScore,score);
					this.PlaceOnRandomCover(i);
					// console.log(coveredScore,score);

					this.state[i] = (this.state[i] == 1) ? -1 : 1;
					this.onExtrema = true;
				}
			} else {
				this.onExtrema = false;
			}
			this.pastParams[i] = {};
			for (let param in this.params[i]) {this.pastParams[i][param] = this.params[i][param];}
			// this.oldScore = score;
		}
	}
	FindBestCover(withApply = false) {
		let best = -Infinity, x = [], ind = -1;
		for (let i = 0; i < this.cover.length; i++) {
			if (this.cover[i][2] > best) {
				best = this.cover[i][2];
				x = this.cover[i][0];
				ind = i;
				// console.log(i);
			}
		}
		// this.cover.sort((a,b)=>{return a[2]-b[2]});
		// if (withApply) {
		// 	let ind = 0;
		// 	for (let param in this.params) {
		// 		this.params[param] = x[ind];
		// 		ind++
		// 	}
		// }
		return [x,best,ind];
	}
}
function sigma(x) {
    return 1 / (1 + Math.pow(Math.E, -x));
}
function tangent(x) {
    return (Math.pow(Math.E, 2 * x) - 1) / (Math.pow(Math.E, 2 * x) + 1);
}
function relu(x) {
    if (x < 0) return 0;
    return x;
}
function reluDer(x) {
    return (x>=0) ? 1:0;
}

function reluLeak(x) {
    if (x < 0) return x/4;
    return x;
}
function reluLeakDer(x) {
    return (x>=0) ? 1:0.25;
}

function EqToArray(func, len) {
    let arr = new Array(len).fill(0).map((v,i)=>{return func(i/len)});
    return arr;
}

class Input {
    constructor(targetNode, Weight = undefined, randomnessRadius = 1) {
        this.node = targetNode;
        if (Weight != undefined)
            this.weight = Weight;
        else
            this.weight = (-1 + Math.random() * 2) * randomnessRadius;
        
        this.conductivity = 1;
        
        this.linkedInputs = [];
    }
    SetNewLinkedInputs(Inputs) {
        // this.weight = Math.random();
        for (let i = 0; i < Inputs.length; i++) {
            Inputs[i].weight = this.weight;
            this.linkedInputs.push(Inputs[i]);
        }
    }
}
class Node {
    constructor(Value = Math.random(), Inputs = []) {
        this.value = Value;
        this.inputs = Inputs;
        this.addicted = [];
        this.deriative = 0;

        this.posInArray = [];

        this.moment1 = [];
        this.moment2 = [];

        this.isBias = false;

        this.clearValue = 1;

        // this.activationFunc = (x) => {return (x);};
        // this.activationFuncDer = (x) => {return 1;}; 

        this.activationFunc = (x) => {return sigma(x);};
        this.activationFuncDer = (x) => {return x*(1-x);}; // this.activationFuncDer = (x) => {return Math.E**-x/((1+Math.E**-x)**2);};
        
        // this.activationFunc = (x) => {return relu(x);};
        // this.activationFuncDer = (x) => {return reluDer(x);};
    }
    addInput(node, randomness = 1, weight = undefined) {
        let inp = new Input(node, weight, randomness);
        this.inputs.push(inp);
        node.addicted.push(new Input(this, inp.weight, randomness));
        return inp;
    }
    getOutput() {
        
    }
}
class NeuroWeb {
    constructor(Nodes = [], ExpectedOutputs = []) { // 1. Nodes = Array of Node class 2. Nodes is 2 dimensional array
        this.nodes = Nodes;
        for (let i = 0; i < this.nodes.length; i++) {
            for (let j = 0; j < this.nodes[i].length; j++) {
                this.nodes[i][j].posInArray = [i, j];
            }
        }
        this.expectedOutputs = ExpectedOutputs;
        this.outputLayerIndex = this.nodes.length - 1;

        this.isRecurrent = false;
        this.recurrentNodesCount = 0;

        this.activationFunc = "";

        this.getOutputs();

        this.asyncTrainerId = null;

        this._debugParameters = {
            Position: [0,0],
            NodesRadius: 30,
            MaxYNodesCount: 5,
            MaxXNodesCount: 8,
            _MaxPosition: []
        };
        document.addEventListener("mousedown", (e)=>{
            this._OnMouseClick(e.offsetX,e.offsetY);
        });
    }
    getOutputs(inps = -1) {
        this.setInputs(inps);
        for (let i = 0; i < this.nodes.length; i++) {
            for (let j = 0; j < this.nodes[i].length; j++) {
                if (this.nodes[i][j].isBias) {
                    this.nodes[i][j].value = 1;
                    continue;
                }

                if (this.nodes[i][j].inputs.length > 0) {
                    // this.nodes[i][j].clearValue = 0;
                    this.nodes[i][j].value = 0;
                } else continue;
                for (let k = 0; k < this.nodes[i][j].inputs.length; k++) {
                    const inp = this.nodes[i][j].inputs[k];
                    this.nodes[i][j].value += (inp.node.value * inp.weight);
                    // this.nodes[i][j].clearValue +=(inp.node.clearValue * inp.weight)
                }
                this.nodes[i][j].value = this.nodes[i][j].activationFunc(this.nodes[i][j].value);
            }
        }
        return this.nodes[this.nodes.length - 1];
    }
    setInputs(arr) {
        if (typeof(arr) == "function") {
            let biasCount = 0;
            for (let i = 0; i < this.nodes[0].length; i++) {
                if (this.nodes[0][i].isBias)
                    biasCount++;
            }
            for (let i = 0; i < this.nodes[0].length - biasCount; i++) {
                let t = i / (this.nodes[0].length - 1 - biasCount);
                this.nodes[0][i].value = arr(t);
            }
            return;
        }
        if (Array.isArray(arr)) {
            for (let i = 0; i < this.nodes[0].length; i++) {
                if (arr.length/this.nodes[0].length < 1) {
                    if (i < arr.length) {
                        this.nodes[0][i].value = arr[i];
                    } else {
                        this.nodes[0][i].value = 0;
                    }
                } else {
                    let ind = Math.floor(i*(arr.length/this.nodes[0].length));
                    this.nodes[0][i].value = arr[ind];
                }
            }
        }
    }
    setInput(i, v) {
        this.nodes[0][i].value = v;
        this.getOutputs();
    }
    getSumError(expectedOuts = this.expectedOutputs) {
        let sum = 0;
        for (let i = 0; i < this.nodes[this.nodes.length - 1].length; i++) {
            sum += Math.pow(expectedOuts[i] - this.nodes[this.nodes.length - 1][i].value, 2);
        }
        return Math.sqrt(sum);
    }
    getWeightsDeriatives(initDeriatives = []) {
        // for (let row = this.nodes.length-1; row >= 0; row--) {
        //     for (let n = 0; n < this.nodes[row].length; n++) {
        //         this.nodes[row][n].deriative = [];
        //         for (let o = 0; o < this.nodes[this.nodes.length-1].length; o++) {
        //             this.nodes[row][n].deriative.push(1);
        //             // this.nodes[row][n].deriative.push((this.expectedOutputs[o]-this.nodes[this.nodes.length-1][o].value));
        //         }
        //     }
        // }
        
        // /*if (typeof InitDeriatives !== "undefined") {
            
        //     for (let out = 0; out < this.nodes[this.nodes.length-1].length; out++)
        //         for (let o = 0; o < InitDeriatives.length; o++) 
        //                 this.nodes[this.nodes.length-1][out].deriative[o] = (InitDeriatives[out][o]);
                
        // }*/
        
        // for (let out = 0; out < this.nodes[this.nodes.length-1].length; out++) {
        //     for (let i = 0; i < this.nodes[this.nodes.length-1][out].inputs.length; i++) {
        //         // if (this.activationFunc == "sigma")
        //             // this.nodes[this.nodes.length-1][out].inputs[i].node.deriative[out] = this.nodes[this.nodes.length-1][out].inputs[i].weight*(this.nodes[this.nodes.length-1][out].value*(1-this.nodes[this.nodes.length-1][out].value))*(this.nodes[this.nodes.length-1][out].inputs[i].node.value*(1-this.nodes[this.nodes.length-1][out].inputs[i].node.value));
        //         // else if (this.activationFunc == "relu") {
        //         //     this.nodes[this.nodes.length-1][out].inputs[i].node.deriative[out] = this.nodes[this.nodes.length-1][out].inputs[i].weight;
        //         //     if (this.nodes[this.nodes.length-1][out].inputs[i].weight < 0) this.nodes[this.nodes.length-1][out].inputs[i].node.deriative[out] = 0;
        //         // } else
        //         //     this.nodes[this.nodes.length-1][out].inputs[i].node.deriative[out] = this.nodes[this.nodes.length-1][out].inputs[i].weight;
        //         this.nodes[this.nodes.length-1][out].inputs[i].node.deriative[out] = this.nodes[this.nodes.length-1][out].inputs[i].weight * this.nodes[this.nodes.length-1][out].activationFuncDer(this.nodes[this.nodes.length-1][out].inputs[i].node.value)//*this.nodes[this.nodes.length-1][out].inputs[i].node.activationFuncDer(this.nodes[this.nodes.length-1][out].inputs[i].node.value);
        //         // this.nodes[this.nodes.length-1][out].inputs[i].node.deriative[out] = this.nodes[this.nodes.length-1][out].inputs[i].weight;
        //         this.nodes[this.nodes.length-1][out].inputs[i].node.newDeriative = this.nodes[this.nodes.length-1][out].newDeriative * this.nodes[this.nodes.length-1][out].inputs[i].weight * this.nodes[this.nodes.length-1][out].activationFuncDer(this.nodes[this.nodes.length-1][out].inputs[i].node.value);
        //         // console.log(outputsDelta[out]);
                
        //     }
        //     for (let row = this.nodes.length-1; row >= 0; row--) {
                
        //         for (let n = 0; n < this.nodes[row].length; n++) {
        //             // this.nodes[row][n].deriative[out] *= 2 * this.nodes[row][n].clearValue;
        //             // if (this.activationFunc == "sigma")
        //                 // this.nodes[row][n].deriative[out] *= this.nodes[row][n].value*(1-this.nodes[row][n].value);
        //             // else if (this.activationFunc == "relu") {
        //             //     if (this.nodes[row][n].value < 0) this.nodes[row][n].deriative[out] = 0;
        //             // }                    
        //             this.nodes[row][n].deriative[out] *= this.nodes[row][n].activationFuncDer(this.nodes[row][n].value);
        //             //this.nodes[row][n].newDeriative *= this.nodes[row][n].activationFuncDer(this.nodes[row][n].value);
        //             for (let i = 0; i < this.nodes[row][n].inputs.length; i++) {
        //                 this.nodes[row][n].inputs[i].node.deriative[out] += this.nodes[row][n].deriative[out] * this.nodes[row][n].inputs[i].weight;
        //                 //this.nodes[row][n].inputs[i].node.newDeriative += this.nodes[row][n].newDeriative * this.nodes[row][n].inputs[i].weight;
        //             }
        //         }
        //     }
        // }
        
        for (let n = 0; n < this.nodes[this.nodes.length-1].length; n++) {
            this.nodes[this.nodes.length-1][n].deriative = initDeriatives[n];
        }
        for (let row = this.nodes.length-1; row >= 0; row--) {
            
            for (let n = 0; n < this.nodes[row].length; n++) {        
                if (row < this.nodes.length-1)
                    this.nodes[row][n].deriative *= this.nodes[row][n].activationFuncDer(this.nodes[row][n].value); 
                for (let i = 0; i < this.nodes[row][n].inputs.length; i++) {
                    this.nodes[row][n].inputs[i].node.deriative = this.nodes[row][n].deriative * this.nodes[row][n].inputs[i].weight;
                }
            }
        }
    }
    normalizeWeights(coff = 1) {
        for (let i = 0; i < this.nodes.length; i++) {
            for (let j = 0; j < this.nodes[i].length; j++) {
                let L = 0;
                for (let k = 0; k < this.nodes[i][j].inputs.length; k++) {
                    L += (this.nodes[i][j].inputs[k].weight)**2;
                }
                L=(L**0.5);
                
                if (L > coff)
                    for (let k = 0; k < this.nodes[i][j].inputs.length; k++) {
                        this.nodes[i][j].inputs[k].weight /= L;
                    }
            }
        }
    }
    GetOutsSubtraction(expected) {return this.nodes[this.nodes.length-1].map((v,i)=>(expected[i]-v.value));}
    GradientStep(deltaError, initDeriatives = null) {
        if (initDeriatives == null) initDeriatives = this.expectedOutputs;
        else this.expectedOutputs = initDeriatives;

        this.getWeightsDeriatives(initDeriatives);
        
        for (let i = 0; i < this.nodes.length; i++) {
            for (let j = 0; j < this.nodes[i].length; j++) {
                for (let k = 0; k < this.nodes[i][j].inputs.length; k++) {

                    let v = this.nodes[i][j].deriative*this.nodes[i][j].inputs[k].node.value * deltaError;; //  * this.nodes[i][j].inputs[k].conductivity //   / (0.5+E)
                    if (isNaN(v)) return;
                    
                    if (this.nodes[i][j].inputs[k].linkedInputs.length > 0)
                        v/=this.nodes[i][j].inputs[k].linkedInputs.length;

                    this.nodes[i][j].inputs[k].weight += v;
                    for (let l = 0; l < this.nodes[i][j].inputs[k].linkedInputs.length; l++) {
                        let linkInps = this.nodes[i][j].inputs[k].linkedInputs[l];

                        linkInps.weight += v;
                    }

                    // let cond = Math.pow(2,-conductivityDecrease*(v**2));
                    // this.nodes[i][j].inputs[k].conductivity *= cond;

                    // // this.nodes[i][j].inputs[k].weight += F/(2*E) * deltaError * this.nodes[i][j].inputs[k].node.value* this.nodes[i][j].inputs[k].conductivity;
                    // let v = F*this.nodes[i][j].inputs[k].node.value;
                    // this.nodes[i][j].inputs[k].weight += v * deltaError * this.nodes[i][j].inputs[k].conductivity;
                    // for (let l = 0; l < this.nodes[i][j].inputs[k].linkedInputs.length; l++) {
                    //     let linkInps = this.nodes[i][j].inputs[k].linkedInputs[l];

                    //     linkInps.weight += v * deltaError;
                    // }

                }
            }
        }
    }
    Adam(Iterations, LearningRate, Targets, beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8, optionalDeriatives = null) {
        let m1 = [];
        let m2 = [];
        for (let i = 0; i < this.nodes.length; i++) {
            m1.push([]);
            m2.push([]);
            for (let j = 0; j < this.nodes[i].length; j++) {
                m1[i].push([]);
                m2[i].push([]);
                for (let k = 0; k < this.nodes[i][j].inputs.length; k++) {
                    m1[i][j].push(0);
                    m2[i][j].push(0);
                }
            }
        }
        
        let t = 0;
        while (t < Iterations) {
            t++;
            for (let i = 0; i < Targets.length; i++) {
                this.getOutputs(Targets[i][0]);
                if (!optionalDeriatives) this.getWeightsDeriatives(this.GetOutsSubtraction(Targets[i][1]));
                else this.getWeightsDeriatives(optionalDeriatives);
                for (let i = 0; i < this.nodes.length; i++) {
                    for (let j = 0; j < this.nodes[i].length; j++) {
                        for (let k = 0; k < this.nodes[i][j].inputs.length; k++) {
                            let g = this.nodes[i][j].deriative*this.nodes[i][j].inputs[k].node.value;
                            m1[i][j][k] = beta1*m1[i][j][k] + (1-beta1)*g;
                            m2[i][j][k] = beta2*m2[i][j][k] + (1-beta2)*g**2;
                            let mHat = m1[i][j][k]/(1-beta1**t);
                            let vHat = m2[i][j][k]/(1-beta2**t);
                            this.nodes[i][j].inputs[k].weight += LearningRate*mHat/(vHat**0.5+epsilon);
                        }
                    }
                }
            }
            //console.log(this.getSumError(Targets[0][1]));
            
        }

    }
    Reset(Random = 1) {
        for (let i = 0; i < this.nodes.length; i++) {
            for (let j = 0; j < this.nodes[i].length; j++) {
                for (let k = 0; k < this.nodes[i][j].inputs.length; k++) {
                    this.nodes[i][j].inputs[k].weight = (-1 + Math.random() * 2)*Random;
                }
            }
        }
        this.getOutputs();
    }
    ResetConductivity() {
        for (let i = 0; i < this.nodes.length; i++) {
            // E = this.getSumError(expectedOutputs);
            for (let j = 0; j < this.nodes[i].length; j++) {
                for (let k = 0; k < this.nodes[i][j].inputs.length; k++) {
                    this.nodes[i][j].inputs[k].conductivity = 1;
                }
            }
        }
    }
    copyWithOffset(DistortionRadius = 0) {
        let nodes = [];

        for (let i = 0; i < this.nodes.length; i++) {
            nodes.push([]);
            for (let j = 0; j < this.nodes[i].length; j++) {
                let N = new Node(this.nodes[i][j].value);
                if (this.nodes[i][j].isBias) N.isBias = true;
                N.activationFunc = this.nodes[i][j].activationFunc;
                N.activationFuncDer = this.nodes[i][j].activationFuncDer;
                nodes[i].push(N);
            }
        }
        for (let i = 0; i < this.nodes.length; i++) {
            for (let j = 0; j < this.nodes[i].length; j++) {
                let inps = [];
                for (let k = 0; k < this.nodes[i][j].inputs.length; k++) {
                    let identifier = this.nodes[i][j].inputs[k].node.posInArray;
                    let w = (this.nodes[i][j].inputs[k].weight + (-1 + Math.random() * 2) * DistortionRadius);
                    inps.push(new Input(nodes[identifier[0]][identifier[1]], w)); //(this.nodes[i][j].inputs[k].weight + (-1 + Math.random() * 2) * DistortionRadius)
                    // console.log(sigma(this.nodes[i][j].inputs[k].weight + (-1 + Math.random() * 2) * DistortionRadius));
                }
                nodes[i][j].inputs = inps;
            }
        }
        return new NeuroWeb(nodes);
    }
    translateFromData(data, step = 1) {
        if (data.length < this.nodes[0].length) {
            for (let i = 0; i < this.nodes[0].length; i++) {
                if (!this.nodes[0][i].isBias)
                    this.nodes[0][i].value = data[i].value;
            }
            this.getOutputs();
        } else {
            for (let i = 0; i < this.nodes[0].length && i < data.length; i += step) {
                this.nodes[0][i].value = data[i].value;
            }
        }
    }
    translateFromEquation(equation, minX = 0, maxX = 1) {
        let biasCount = 0;
        for (let i = 0; i < this.nodes[0].length; i++) {
            if (this.nodes[0][i].isBias)
                biasCount++;
        }
        for (let i = 0; i < this.nodes[0].length - biasCount; i++) {
            let t = i / (this.nodes[0].length - 1 - biasCount);
            this.nodes[0][i].value = equation(minX + t * (maxX - minX));
        }
    }
    _GetPositionByIndexes(i, j, r = 20) {
        return [i * (r + 20) * 2, j * (r + 5) * 2];
    }
    _OnMouseClick(x,y) {
        for (let i = 0; i < this.nodes[0].length; i++) {
            let nPos = this._GetPositionByIndexes(0,i,20);
            nPos = [nPos[0]+this._debugParameters.Position[0],nPos[1]+this._debugParameters.Position[1]]; 
            
            let delta = [x-nPos[0],y-nPos[1]];
            if (math.magnitude(delta) < this._debugParameters.NodesRadius) {
                // this.nodes[0][i].value = Number(prompt("Set value"));
                this.getOutputs();
            }
        }
    }
    _DebugDraw(StartPosition = { x: 0, y: 0 }, r = 20) {
        // this.getOutputs();
        this._debugParameters.Position = [StartPosition.x,StartPosition.y]
        this._debugParameters.NodesRadius = r;

        let maxNodesCount = this._debugParameters.MaxYNodesCount;
        let maxNodesCountX = this._debugParameters.MaxXNodesCount;

        let lastI = 0;
        let maxY = 0;
        for (let i = 0; i < this.nodes.length; i++) {
            let indexOffset = 0;
            lastI = i;
            if (i > maxNodesCountX) {
                lastI = maxNodesCountX+2;
                i = this.nodes.length-1;

                for (let j = 0; j < this.nodes[this.nodes.length-1].length; j++) {
                    let pos = this._GetPositionByIndexes(lastI,j-indexOffset,r);
                    let v = this.nodes[this.nodes.length-1][j].value*255;
                    let clr = inRgb(v,-v,0,1);
                    d.circle(StartPosition.x+pos[0],StartPosition.y+pos[1],r,clr,clr,1);

                    let hiddenPos = this._GetPositionByIndexes(maxNodesCountX+1,j,r);
                    d.rect(StartPosition.x+hiddenPos[0],StartPosition.y+hiddenPos[1],2,2,"white","white");
                    d.rect(StartPosition.x+hiddenPos[0]-15,StartPosition.y+hiddenPos[1],2,2,"white","white");
                    d.rect(StartPosition.x+hiddenPos[0]+15,StartPosition.y+hiddenPos[1],2,2,"white","white");
                }
                continue;
            }
            d.txt(this.nodes[i].length.toString(), StartPosition.x+this._GetPositionByIndexes(lastI,0,r)[0]-8, StartPosition.y-r,"12px Arial","white");
            for (let j = 0; j < this.nodes[i].length; j++) {
                // Nodes
                if (j > maxNodesCount) {
                    indexOffset = this.nodes[i].length-1-maxNodesCount-2;
                    j = this.nodes[i].length-1;

                    let hiddenPos = this._GetPositionByIndexes(i,maxNodesCount+1,r);
                    d.rect(StartPosition.x+hiddenPos[0],StartPosition.y+hiddenPos[1],2,2,"white","white");
                    d.rect(StartPosition.x+hiddenPos[0],StartPosition.y+hiddenPos[1]-10,2,2,"white","white");
                    d.rect(StartPosition.x+hiddenPos[0],StartPosition.y+hiddenPos[1]+10,2,2,"white","white");
                }

                let pos = this._GetPositionByIndexes(lastI,j-indexOffset,r);
                let v = this.nodes[lastI][j].value*255;
                let clr = inRgb(v,-v,0,1);
                d.circle(StartPosition.x+pos[0],StartPosition.y+pos[1],r,clr,clr,1);
                d.circle(StartPosition.x+pos[0],StartPosition.y+pos[1],r,clr,clr,1);
                if (pos[1] > maxY) maxY = pos[1];

                // Lines
                for (let k = 0; k < this.nodes[lastI][j].inputs.length; k++) {
                    const line = this.nodes[lastI][j].inputs[k];

                    let v = line.weight*255*5;
                    let clr = inRgb(v,-v,0,0.7);
                 
                    if (line.node.posInArray[1] > maxNodesCount) {
                        if (line.node.posInArray[1] == this.nodes[line.node.posInArray[0]].length-1) {
                            let targetLinePos = this._GetPositionByIndexes(line.node.posInArray[0],maxNodesCount+2,r);

                            d.line(StartPosition.x+pos[0],StartPosition.y+pos[1], StartPosition.x+targetLinePos[0],StartPosition.y+targetLinePos[1], clr);
                        }
                        continue;
                    }
                    let targetLinePos = this._GetPositionByIndexes(line.node.posInArray[0],line.node.posInArray[1],r);

                    let dst = math.distanceToLine([input.mouse.x,input.mouse.y],[StartPosition.x+pos[0],StartPosition.y+pos[1]],[StartPosition.x+targetLinePos[0],StartPosition.y+targetLinePos[1]]);
                    if (dst[0] < 5) {
                        d.txt(line.weight.toFixed(3),input.mouse.x,input.mouse.y,"16px Arial","white");
                    }

                    d.line(StartPosition.x+pos[0],StartPosition.y+pos[1], StartPosition.x+targetLinePos[0],StartPosition.y+targetLinePos[1], clr);
                }
            }
        }

        let rightWebBound = this._GetPositionByIndexes(lastI+1,0)[0];
        let rightBlockPos = [
            StartPosition.x+rightWebBound-r*2,StartPosition.y-r
        ];
        d.line(StartPosition.x+rightWebBound-r*2,StartPosition.y-r, StartPosition.x+rightWebBound-r*2, StartPosition.y+maxY+r,"white",1);

        this._debugParameters._MaxPosition = [StartPosition.x+rightWebBound-r*2,StartPosition.y+maxY+r];

        
        // d.txt("Inputs",rightBlockPos[0]+5,rightBlockPos[1]+16,"16px Arial","white");
        // for (let i = 0; i < this.nodes[0].length; i++) {
        //     d.txt(this.nodes[0][i].value.toFixed(3),rightBlockPos[0]+5+8,rightBlockPos[1]+16+(i+1)*16,"14px Arial","white");
        // }
        d.txt("Outputs",rightBlockPos[0]+5,rightBlockPos[1]+16,"16px Arial","white");
        for (let i = 0; i < this.nodes[this.nodes.length-1].length; i++) {
            d.txt(this.nodes[this.nodes.length-1][i].value.toFixed(3),rightBlockPos[0]+5+8,rightBlockPos[1]+16+(i+1)*16,"14px Arial","white");
        }
        d.txt("Error: " + this.getSumError().toFixed(3),rightBlockPos[0]+5,rightBlockPos[1]+16+((this.nodes[this.nodes.length-1].length+1)*16)+16,"16px Arial","white");
        d.txt("Expected values: ",rightBlockPos[0]+5,rightBlockPos[1]+16+((this.nodes[this.nodes.length-1].length+1)*16)+16+16,"16px Arial","white");
        for (let i = 0; i < this.expectedOutputs.length; i++) {
            d.txt(""+this.expectedOutputs[i].toFixed(3),rightBlockPos[0]+5+8,rightBlockPos[1]+16+((this.nodes[this.nodes.length-1].length+1)*16)+16+16+16+16*i,"16px Arial","white");
        }
    }
    _DebugDrawOutputs(StartPosition = { x: 0, y: 0 }, DrawInps = true) {
        this.getOutputs();
        let maxLen = 0;
        ctx.font = "12px Arial";
        let maxInpLen = 0;
        if (DrawInps) {
            d.txt("Inputs> ", StartPosition.x,StartPosition.y,"","white");
            for (let i = 0; i < this.nodes[0].length; i++) {
                let len = ctx.measureText(this.nodes[0][i].value.toFixed(3)).width;
                d.txt(this.nodes[0][i].value.toFixed(3),StartPosition.x,StartPosition.y+i*12+12,"","White");
                if (maxInpLen < len) maxInpLen = len;
            }
        }
        d.txt("Outputs> ", StartPosition.x+50,StartPosition.y,"","white");
        for (let i = 0; i < this.nodes[this.nodes.length - 1].length; i++) {
            d.txt(this.nodes[this.nodes.length - 1][i].value.toFixed(3), StartPosition.x+maxInpLen+20, StartPosition.y + 12 * (i + 1), "", "white");
            let len = ctx.measureText(toString(this.nodes[this.nodes.length - 1][i].value.toFixed(3))).width;
            if (len > maxLen) {
                maxLen = len;
            }
        }
        ctx.font = "";
        d.line(StartPosition.x + maxLen + maxInpLen, StartPosition.y, StartPosition.x + maxLen + maxInpLen, StartPosition.y + 12 * this.nodes[this.nodes.length - 1].length, "white", 1);
        d.txt(this.getSumError(this.expectedOutputs).toFixed(3), StartPosition.x + maxLen + maxInpLen + 5, StartPosition.y + 6 * this.nodes[this.nodes.length - 1].length + 6);
    }
}
class WebTrainer {
    #graphDataIndices = [];
    constructor(TargetWeb, DebugGraph = null) {
        this.web = TargetWeb;

        this.inputData = []; // [[1,2,3,4],[4,3,2,1]] - Example
        this.outputData = []; // [[1,0], [0,1]] - Example
        this.deltaError = 0.01;
        
        this.bestConfiguration = {
            minError: Infinity,
            inputsConfiguration: []
        };

        this.graph = DebugGraph;
        this.errorGraphInd = this.CreateNewDataInGraph("Training progress","rgba(0,255,0,0.8)");
    }
    AddTrainData(inputs, outputs = inputs) {
        if (typeof(inputs) == "function") {
            let biasCount = 0;
            for (let i = 0; i < this.web.nodes[0].length; i++) {
                if (this.web.nodes[0][i].isBias)
                    biasCount++;
            }
            this.inputData.push([]);
            for (let i = 0; i < this.web.nodes[0].length - biasCount; i++) {
                let t = i / (this.web.nodes[0].length - 1 - biasCount);
                this.inputData[this.inputData.length-1].push(inputs(t));
            }
        } else
            this.inputData.push(inputs);
        if (typeof(outputs) == "function") {
            this.outputData.push([]);
            for (let i = 0; i < this.web.nodes[this.web.nodes.length-1].length; i++) {
                let t = i / (this.web.nodes[this.web.nodes.length-1].length - 1);
                this.outputData[this.outputData.length-1].push(outputs(t));
            }
        } else
            this.outputData.push(outputs);

            // else if (typeof(outputs) == "undefined") {
            //     this.outputData.push(inputs);
            // }
            
        // this.CreateNewDataInGraph("Err"+(this.#graphDataIndices.length-1),inRgb((0.5+Math.random()*0.5)*255,(0.5+Math.random()*0.5)*255,(0.5+Math.random()*0.5)*255));
    }
    TrainWeb(deltaError = null, epochs = 1) {
        if (deltaError == null) deltaError = this.deltaError;
        for (let T = 0; T < epochs; T++) {
            let error = 0;
            for (let i = 0; i < this.inputData.length; i++) {
                this.web.setInputs(this.inputData[i]);//NewGraph.ParseEquationAsData((x)=>{return Math.random()},this.inputData[i].length,[0,1])
                // this.web.setInputs(NewGraph.ParseEquationAsData((x)=>{return Math.random()},this.inputData[i].length,[0,1]));//
                this.web.oneIterTrain(deltaError,this.outputData[i]);

                let E = this.web.getSumError();
                error += E;
                
                //if (this.graph != null)
                //    this.AddDataToGraph(E,this.errorGraphInd);
            }
            error/=this.inputData.length;
            // console.log(error);

            if (this.graph != null)
                this.AddDataToGraph(error,0);
        }
    }
    GetErrorsData() {
        let startValues = [];
        for (let i = 0; i < this.web.nodes[0].length; i++) {
            startValues.push(this.web.nodes[0][i].value);
        }

        let data = [];
        for (let i = 0; i < this.inputData.length; i++) {
            this.web.getOutputs(this.inputData[i]);
            data.push(this.web.getSumError(this.outputData[i]));
        }

        this.web.getOutputs(startValues);
        return data;
    }
    GetOutputsData() {
        let startValues = [];
        for (let i = 0; i < this.web.nodes[0].length; i++) {
            startValues.push(this.web.nodes[0][i].value);
        }

        let data = [];
        for (let i = 0; i < this.inputData.length; i++) {
            let nwData = this.web.getOutputs(this.inputData[i]);
            data.push([]);
            for (let j = 0; j < nwData.length; j++) {
                data[data.length-1].push(nwData[j].value);
            }
        }

        this.web.getOutputs(startValues);
        return data;
    }
    AddDataToGraph(Data, Index, XOffset = 0) {
        if (Index < 0 || Index > this.#graphDataIndices.length-1) return null;
        this.graph.data[this.#graphDataIndices[Index]].data.push({x:(this.graph.data[this.#graphDataIndices[Index]].data.length-1)/1000 + XOffset,y:Data});
        // if (this.graph.data[this.#graphDataIndices[Index]].data.length > 1000) {
        //     this.graph.data[this.#graphDataIndices[Index]].data.splice(0,1);
        //     for (let i = 0; i < this.graph.data[this.#graphDataIndices[Index]].data.length; i++) {
        //         this.graph.data[this.#graphDataIndices[Index]].data[i].x -= 1/1000;
        //     }
        // }
    }
    CreateNewDataInGraph(Name, Color) {
        if (this.graph != null) {
            this.graph.data[this.graph.data.length] = {name:Name,color:Color,data:[]};
            this.#graphDataIndices[this.#graphDataIndices.length] = this.graph.data.length-1;
        } else {return null;}

        return this.#graphDataIndices[this.#graphDataIndices.length-1];
    }
}
function createNWColumnOptions(data) {
    let defaultOptions = {
        bias:0,
        connectType:"full",
        convolutionRadius:1,
        activationFunc: (x)=>{return (tangent(x));},
        activationFuncDer: (x)=>{return 1-x**2;},
        // activationFunc: (x)=>{return (sigma(x));},
        // activationFuncDer: (x)=>{return x*(1-x);},
        // activationFunc: (x)=>{return relu(x);},
        // activationFuncDer: (x)=>{return reluDer(x);},//x*(1-x)
        // activationFunc: (x)=>{return (x);},
        // activationFuncDer: (x)=>{return (1);},//x*(1-x)
        randomness: 1,
    };
    if (typeof (data) == "undefined")
        data = defaultOptions;
    for (option in defaultOptions) {
        if (typeof (data[option]) == "undefined") {
            data[option] = defaultOptions[option];
        }
    }
    return data;
}
function GenerateNW(columnsData, options = []) {
    options = options.flat(Infinity);
    for (let i =  0; i < columnsData.length; i++) {
        options[i] = createNWColumnOptions(options[i]);
    }

    let nodes = [];
    for (let i = 0; i < columnsData.length; i++) {
        let column = [];
        let isConvoluted = false;
        if (options[i].connectType == "convolution") {
            if (columnsData[i] == columnsData[i-1]) {
                isConvoluted = true;
            } else {
                options[i] = createNWColumnOptions({
                    connectType: "compress",
                    activationFunc: (x)=>{return x;},
                    activationFuncDer: (x)=>{return 0;},
                    bias:options[i].bias,
                    randomness:options[i].randomness
                });
                columnsData.splice(i,0,columnsData[i]);
                options.splice(i+1,0,createNWColumnOptions({
                    connectType: "convolution",
                    activationFunc: (x)=>{return relu(x);},
                    activationFuncDer: (x)=>{return reluDer(x);},
                    randomness:options[i].randomness
                }));
            }
        }
        for (let j = 0; j < columnsData[i]+options[i].bias; j++) {
            let N = new Node();
            if (j >= columnsData[i]) {
                N.value = 1;
                N.isBias = true;
            }
            N.activationFunc = options[i].activationFunc;
            N.activationFuncDer = options[i].activationFuncDer;
            if (i == columnsData.length-1) {
                // N.activationFunc = (x) => {return x;};
                // N.activationFuncDer = (x) => {return 1;};
            }

            if (options[i].connectType == "full" && j < columnsData[i]) {
                for (let k = 0; i > 0 && k < nodes[i-1].length; k++) {
                    N.addInput(nodes[i-1][k],options[i-1].randomness);
                }
            } else if (options[i].connectType == "convolution" && j < columnsData[i]) {
                for (let r = -1*options[i].convolutionRadius; r <= 1*options[i].convolutionRadius; r++) {
                    if (j+r < 0 || j+r >= nodes[i-1].length) continue;
                    N.addInput(nodes[i-1][j+r],options[i-1].randomness);
                }
                if (options[i-1].bias > 0) {
                    N.addInput(nodes[i-1][nodes[i-1].length-1],options[i].randomness);
                }
            } else if (options[i].connectType == "compress" && j < columnsData[i]) {
                let step = Math.ceil(columnsData[i-1]/columnsData[i]);
                isConvoluted = true;
                for (let r = 0; r < step; r++) {
                    if (r+j*step >= nodes[i-1].length) continue;
                    let inp = N.addInput(nodes[i-1][r+j*step],options[i].randomness);
                    inp.weight = columnsData[i]/columnsData[i-1];
                    inp.conductivity = 0;
                }
            }
            column.push(N);
        }

        nodes.push(column);
        
        if (isConvoluted) {
            /*for (let j = options[i].convolutionRadius; j < columnsData[i]-options[i].convolutionRadius; j++) {
                for (let k = 0; k <= options[i].convolutionRadius*2; k++) {
                    let links = [];
                    for (let l = options[i].convolutionRadius; l < columnsData[i]-options[i].convolutionRadius; l++) {
                        if (j==l) continue;
                        if (typeof (nodes[i][l].inputs[k]) == "undefined") continue;
                        links.push(nodes[i][l].inputs[k]);
                    }
                    if (typeof (nodes[i][j].inputs[k]) == "undefined") continue;
                    nodes[i][j].inputs[k].SetNewLinkedInputs(links);
                }
            }*/
            for (let j = 0; j < nodes[i].length; j++) {
                for (let k = 0; k < nodes[i][j].inputs.length && !nodes[i][j].isBias; k++) {
                    let links = [];
                    for (let l = 0; l < nodes[i].length; l++) {
                        if (j==l) continue;
                        if (typeof (nodes[i][l].inputs[k]) == "undefined") continue;
                        links.push(nodes[i][l].inputs[k]);
                    }
                    nodes[i][j].inputs[k].SetNewLinkedInputs(links);
                }
            }
        }
    }
    let nw = new NeuroWeb(nodes);
    return nw;
}
function GenerateConvolutionalNeuroWeb(InputsCount, ConvLayersCount, HiddenX, HiddenY, OutputsCount, ConvolutionRadius = 1, DiscretizeRate = 2, Randomness = 1, HaveBias = 0) {
    let nodes = [];
    
    nodes.push([]);
    for (let i = 0; i < InputsCount; i++) {
        nodes[0].push(new Node());
    }

    for (let x = 0; x < ConvLayersCount; x++) {
        nodes.push([]);
        let layerNodesCount = InputsCount;
        if (DiscretizeRate > 0)
            layerNodesCount /= (DiscretizeRate**x);
        for (let y = 0; y < layerNodesCount; y++) {
            let N = new Node();
            for (let k = -ConvolutionRadius; k <= ConvolutionRadius; k++) {
                let I = (y+k);
                // if (nodes[nodes.length-1-1].length < DiscretizeRate) continue;
                if (I < 0) {
                    N.addInput(nodes[nodes.length-1-1][0],Randomness);
                    continue;
                }
                if (I > layerNodesCount-1) {
                    N.addInput(nodes[nodes.length-1-1][layerNodesCount-1],Randomness);
                    continue;
                }
                N.addInput(nodes[nodes.length-1-1][I],Randomness);
            }
            N.activationFunc = (x)=>{return relu(x);};
            N.activationFuncDer = (x)=>{return reluDer(x);};
            nodes[nodes.length-1].push(N);
        }
        for (let y = 0; y < nodes[nodes.length-1].length; y++) {
            for (let i = 0; i < nodes[nodes.length-1][y].inputs.length; i++) {
                if (nodes[nodes.length-1][y].isBias) continue;
                let newLinks = [];
                for (let l = 0; l < nodes[nodes.length-1].length; l++) {
                    if (l==y) continue;
                    newLinks.push(nodes[nodes.length-1][l].inputs[i]);
                }
                nodes[nodes.length-1][y].inputs[i].SetNewLinkedInputs(newLinks);
            }
        }

        if (DiscretizeRate > 0) {
            nodes.push([]);
            for (let y = 0; y < InputsCount/(DiscretizeRate**(x+1)); y++) {
                let N = new Node();
                
                for (let k = 0; k < DiscretizeRate; k++) {
                    let inp = N.addInput(nodes[nodes.length-1-1][k+y*DiscretizeRate],0,1/DiscretizeRate,Randomness);
                    inp.conductivity = 0;
                }

                N.activationFunc = (x)=>{return relu(x);};
                N.activationFuncDer = (x)=>{return reluDer(x);};
                N.activationFunc = (x)=>{return (x);};
                N.activationFuncDer = (x)=>{return 0;};
                nodes[nodes.length-1].push(N);
            }
            for (let y = 0; y < nodes[nodes.length-1].length*1; y++) {
                for (let i = 0; i < DiscretizeRate; i++) {
                    let newLinks = [];
                    for (let l = 0; l < nodes[nodes.length-1].length; l++) {
                        if (l==y) continue;
                        newLinks.push(nodes[nodes.length-1][l].inputs[i]);
                    }
                    nodes[nodes.length-1][y].inputs[i].SetNewLinkedInputs(newLinks);
                }
            }
        }
        if (HaveBias) {
            let N = new Node();
            N.isBias = true;
            N.value = 1;
            nodes[nodes.length-1].push(N)
            for (let i = 0; i < nodes[x*2+1].length; i++) {
                nodes[x*2+1][i].addInput(nodes[x*2][nodes[x*2].length-1], Randomness);
            }
        }
    }

    for (let x = 0; x < HiddenX; x++) {
        nodes.push([]);
        for (let y = 0; y < HiddenY; y++) {
            let N = new Node();
            for (let k = 0; k < nodes[nodes.length-1-1].length; k++) {
                N.addInput(nodes[nodes.length-1-1][k],Randomness);
            }
            N.activationFunc = (x) => {return sigma(x);};
            N.activationFuncDer = (x) => {return x*(1-x);};
            // N.activationFunc = (x) => {return relu(x);};
            // N.activationFuncDer = (x) => {return reluDer(x);};
            nodes[nodes.length-1].push(N);
        }
    }

    if (OutputsCount > 0)
        nodes.push([]);
    for (let i = 0; i < OutputsCount; i++) {
        let N = new Node();
        for (let k = 0; k < nodes[nodes.length-1-1].length; k++) {
            N.addInput(nodes[nodes.length-1-1][k],Randomness);
        }
        N.activationFunc = (x) => {return relu(x);};
        N.activationFuncDer = (x) => {return reluDer(x);};
        // N.activationFunc = (x) => {return sigma(x);};
        // N.activationFuncDer = (x) => {return x/(1+x);};
        nodes[nodes.length-1].push(N);
    }
    for (let i = 0; i < nodes[nodes.length-1].length; i++) {
        for (let k = 0; k < nodes[nodes.length-1][i].inputs.length; k++) {
            let links = [];
            for (let l = 0; l < nodes[nodes.length-1][i].inputs.length; l++) {
                if (k == l) continue;
                links.push(nodes[nodes.length-1][i].inputs[l]);
            }
            // nodes[nodes.length-1][i].inputs[k].SetNewLinkedInputs(links);
        }
    }

    let NW = new NeuroWeb(nodes);
    NW.getOutputs();

    return NW;
}
