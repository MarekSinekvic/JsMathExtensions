var d = new (class _ {
	constructor() { }

	circle(x, y, r, fillcolor = "black", linecolor = "black", lineWidth = 1, filled = true, angles = [0, 360], direction = false) {
		let startLineWidth = ctx.lineWidth;
		let startColor = { fill: ctx.fillStyle, stroke: ctx.strokeStyle };

		ctx.lineWidth = lineWidth;
		ctx.strokeStyle = linecolor;
		ctx.fillStyle = fillcolor;
		ctx.beginPath();
		ctx.arc(x, y, r, (angles[0] * Math.PI) / 180, (angles[1] * Math.PI) / 180, direction);
		if (filled)
			ctx.fill();
		else
			ctx.stroke();
		ctx.lineWidth = startLineWidth;
		ctx.fillColor = startColor.fill;
		ctx.strokeColor = startColor.stroke;
	}

	rect(x, y, w, h, fillcolor = "black", linecolor = "black", linewidth = 1, filled = true) {
		ctx.fillStyle = fillcolor;
		ctx.strokeStyle = linecolor;
		ctx.lineWidth = linewidth;
		if (filled) {
			ctx.fillRect(x, y, w, h);
		} else {
			ctx.strokeRect(x, y, w, h);
		}
	}

	line(x, y, x1, y1, clr = "black", lineWidth = 1) {
		ctx.fillStyle = clr;
		ctx.strokeStyle = clr;
		ctx.beginPath();
		ctx.lineWidth = lineWidth;
		ctx.moveTo(x, y);
		ctx.lineTo(x1, y1);
		ctx.stroke();
		ctx.fillStyle = "black";
		ctx.strokeStyle = "black";
		ctx.lineWidth = 1;
	}

	ray(x, y, direction, length = 10, clr = "white", withCtxPath = true, inLengthArrowSize = 0.6, inNormalArrowSize = 0.2, lineWidth = ctx.lineWidth) {
		let startLineWidth = ctx.lineWidth;

		direction = math.normalize(direction);

		ctx.lineWidth = lineWidth;

		let normal = { x: -direction[1], y: direction[0] };

		if (withCtxPath) {
			ctx.strokeStyle = clr;
			ctx.beginPath();
		}
		ctx.moveTo(x, y);

		ctx.lineTo(x + direction[0] * length, y + direction[1] * length);

		ctx.lineTo(x + direction[0] * (length * inLengthArrowSize) + normal.x * (inNormalArrowSize * length), y + direction[1] * (length * inLengthArrowSize) + normal.y * (inNormalArrowSize * length));
		ctx.moveTo(x + direction[0] * length, y + direction[1] * length);
		ctx.lineTo(x + direction[0] * (length * inLengthArrowSize) - normal.x * (inNormalArrowSize * length), y + direction[1] * (length * inLengthArrowSize) - normal.y * (inNormalArrowSize * length));

		if (withCtxPath)
			ctx.stroke();

		ctx.lineWidth = startLineWidth;
	}
	dray(x, y, direction, clr = "white", inLengthArrowSize = 0.6, inNormalArrowSize = 0.2, lineWidth = ctx.lineWidth) {
		let startLineWidth = ctx.lineWidth;

		ctx.lineWidth = lineWidth;

		let normal = ({ x: -direction.y, y: direction.x });

		ctx.strokeStyle = clr;
		ctx.beginPath();

		ctx.moveTo(x, y);

		ctx.lineTo(x + direction.x * length, y + direction.y * length);

		ctx.lineTo(x + direction.x * (length * inLengthArrowSize) + normal.x * (inNormalArrowSize * length), y + direction.y * (length * inLengthArrowSize) + normal.y * (inNormalArrowSize * length));
		ctx.moveTo(x + direction.x * length, y + direction.y * length);
		ctx.lineTo(x + direction.x * (length * inLengthArrowSize) - normal.x * (inNormalArrowSize * length), y + direction.y * (length * inLengthArrowSize) - normal.y * (inNormalArrowSize * length));

		ctx.stroke();

		ctx.lineWidth = startLineWidth;
	}

	txt(text, x, y, font, color = "white", isStroke = false) {
		let startFont = ctx.font;
		let startColor = { fill: ctx.fillStyle, stroke: ctx.strokeStyle };

		ctx.font = font;
		ctx.fillStyle = color;
		ctx.strokeStyle = color;

		ctx.fillText(text, x, y);

		ctx.font = startFont;
		ctx.fillColor = startColor.fill;
		ctx.strokeColor = startColor.stroke;
	}
	Matrix(M) {

	}

	clear(clr) {
		canv.style.backgroundColor = clr;
		ctx.clearRect(0, 0, canv.width, canv.height);
	}

	drawImg(img, x, y, w, h, sx = 0, sy = 0, sw = w, sh = h) {
		ctx.drawImage(img, sx, sy, sw, sh, x, y, w, h);
	}
})();
class Slider {
	constructor(Position = [0,0], Label = "test", OnChange = new function (v) { }, SimbolsCountByDot = 2, Length = 120, StartCurrentValue = 50, MinValue = -100, MaxValue = 100, ValueStep = 0) {
		this.label = Label;
		ctx.font = "12px Arial";
		this.labelLength = ctx.measureText(this.label).width;

		this.position = Position;
		this.position[0] += this.labelLength - 4;
		this.length = Length - ctx.measureText(Label).width - ctx.measureText((Math.pow(10, SimbolsCountByDot + 1) - 1).toString()).width;
		this.min = MinValue;
		this.max = MaxValue;
		this.step = ValueStep;

		this.simbolsCountByDot = SimbolsCountByDot;

		this.onChange = OnChange;

		this.stickPosition = StartCurrentValue / MaxValue;

		this.wasClicked = false;

		// document.addEventListener("mousedown", (e)=>{
		// 	this.position[0] + this.length + 5
		// })
		document.addEventListener("mousemove",(e)=>{
			// console.log(e);
			if (e.buttons != 1) return;
			this.#control();
		});
		document.addEventListener("mouseup", (e)=>{
			this.wasClicked = false;
		});
	}
	Draw() {
		ctx.font = "12px Arial";
		d.line(this.position[0], this.position[1], this.position[0] + this.length, this.position[1], "white");
		d.line(this.position[0], this.position[1] - 4, this.position[0], this.position[1] + 4, "white");
		d.line(this.position[0] + this.length, this.position[1] - 4, this.position[0] + this.length, this.position[1] + 4, "white");

		d.circle(this.position[0] + this.stickPosition * this.length, this.position[1], 4, "white", "white", 1, false);
		d.circle(this.position[0] + this.stickPosition * this.length, this.position[1], 3, "black", "black", 1, true);

		d.txt(this.getResult().toFixed(this.simbolsCountByDot), this.position[0] + this.length + 5, this.position[1] + 9 / 3, "9px Arial");
		d.txt(this.label, this.position[0] - this.labelLength - 4, this.position[1] + 5, "", "white");

	}
	#control() {
		let dst = math.distance([input.mouse.x, input.mouse.y], [this.position[0] + this.stickPosition * this.length, this.position[1]]);
		let clickPos = -1;
		if (dst[0] < 4)
			this.wasClicked = true;
		if (Math.abs(dst[2]) < 4 && Math.abs(dst[1]) < this.length) {
			this.wasClicked = true;
			clickPos = (input.mouse.x - this.position[0]) / this.length;
		}
		if (this.wasClicked) {
			if (clickPos != -1) {
				this.stickPosition = clickPos;
				clickPos = -1;
			}
			this.stickPosition = (input.mouse.x - this.position[0]) / this.length;
			if (this.stickPosition > 1) this.stickPosition = 1;
			if (this.stickPosition < 0) this.stickPosition = 0;

			this.onChange(this.getResult());
		}
	}
	setStickPosition(t) {
		this.stickPosition = (this.min - t) / (this.min - this.max);
	}
	isHoverOnSlider(point) {
		let dst = math.distance([point[0], point[1]], [this.position[0] + this.stickPosition * this.length, this.position[1]]);
		if (Math.abs(dst[2]) < 4 && Math.abs(dst[1]) < this.length) {
			return true;
		}
		return false;
	}
	getResult() {
		return (this.min + this.stickPosition * (this.max - this.min));
	}
	SetValue(v) {
		this.stickPosition = v/(this.max-this.min);
	}
}
class Toggle {
	#onClickFuncs;
	constructor(Position, TextSize, StartState = false, TextOnFalse = "FALSE", TextOnTrue = "TRUE") {
		this.position = Position;
		this.size = TextSize;
		this.state = StartState;

		this.textOnFalse = TextOnFalse;
		this.textOnTrue = TextOnTrue;

		this.hotkeyButton = "";

		this.#onClickFuncs = [];

		document.addEventListener("mousedown", (e) => { this.OnClick(e); });
		document.addEventListener("keydown", (e) => {
			if (this.hotkeyButton == "") return;
			if (e.code == this.hotkeyButton) {
				this.SwapState();
			}
		});
	}
	GetBordersRect() {
		let width = 0;
		let startFont = ctx.font;
		ctx.font = this.size + "px Arial";
		if (this.state)
			width = ctx.measureText(this.textOnTrue).width;
		else
			width = ctx.measureText(this.textOnFalse).width;

		ctx.font = startFont;
		return {
			x: this.position[0] - 5,
			y: this.position[1] - 5,
			w: width + 5 * 2,
			h: this.size + 5 * 2
		};
	}
	addOnClickEvent(func) {
		this.#onClickFuncs.push(func);
	}
	SwapState() {
		this.state = !this.state;
		for (let i = 0; i < this.#onClickFuncs.length; i++) {
			this.#onClickFuncs[i](this.state);
		}
	}
	IsMouseHovered() {
		let bordersRect = this.GetBordersRect();
		if (input.mouse.x > bordersRect.x && input.mouse.x < bordersRect.x + bordersRect.w &&
			input.mouse.y > bordersRect.y && input.mouse.y < bordersRect.y + bordersRect.h) {
				return true;
		}
		return false;
	}
	OnClick(e) {
		if (this.IsMouseHovered()) {
			this.SwapState();
		}
	}
	Draw() {
		let bordersRect = this.GetBordersRect();
		d.rect(bordersRect.x, bordersRect.y, bordersRect.w, bordersRect.h, (this.state) ? "rgba(255,255,255,1)" : "black", "white", 1, true);
		d.rect(bordersRect.x, bordersRect.y, bordersRect.w, bordersRect.h, "black", "white", 1, false);

		let startFont = ctx.font;
		ctx.font = this.size + "px Arial";
		let width = 0;
		if (this.state) {
			width = ctx.measureText(this.textOnTrue).width;
			d.txt(this.textOnTrue, this.position[0], this.position[1] + this.size - 3, this.size + "px Arial", "black", false);
		} else {
			width = ctx.measureText(this.textOnFalse).width;
			d.txt(this.textOnFalse, this.position[0], this.position[1] + this.size - 3, this.size + "px Arial", "white", true);
		}
		ctx.font = startFont;
	}
}
class GraphData {
	constructor(Value, XPosition, Color = "rgb(255,255,255)", Label = "") {
		this.value = Value;
		this.xpos = XPosition;
		this.color = Color;
		this.label = Label;
	}
	GetDeriative(x, d = 0.000001) {
		return (this.value(x + d) - this.value(x)) / d;
	}
}
class UIInput {
	#onChangeEvents;
	constructor(Position, Width, FontSize, Value = "") {
		this.position = Position;
		this.width = Width;
		this.height = FontSize+6;
		this.fontSize = FontSize;

		this.value = Value;

		ctx.font = this.fontSize+"px Arial";

		this._isInTarget = false;
		
		this._pixelTextLength = ctx.measureText(this.value).width;
		
		this._textCursorColor = "white";
		this._cursorPosition = 0;
		this._selectionBounds = [-1,-1];
		setInterval(()=>{
			if (this._textCursorColor == "white") {
				this._textCursorColor = "black";
			} else if (this._textCursorColor == "black") {
				this._textCursorColor = "white";
			}
		}, 300);

		this.#onChangeEvents = [];

		document.addEventListener("mousedown", (e)=>{
			let clkPos = [e.offsetX,e.offsetY];
			if ((clkPos[0] > this.position[0] && clkPos[0] < this.position[0]+this.width) &&
				(clkPos[1] > this.position[1] && clkPos[1] < this.position[1]+this.height)) {
				this._isInTarget = true;
			}  else {
				this._isInTarget = false;
			}
			
			if (e.buttons == 2 && this._isInTarget) console.info(this.value);
			if (e.buttons == 1 && this._isInTarget) {
				let localXClkPos = clkPos[0]-this.position[0];
				let minDst = Infinity;
				let ind = -1;

				ctx.font = this.fontSize+"px Arial";
				for (let i = 0; i < this.value.length; i++) {
					let pos = ctx.measureText(this.value.slice(0,i+1)).width;
					let dst = Math.abs(localXClkPos-pos);
					if (dst < minDst) {
						minDst = dst;
						ind = i;
					}
				}

				if (ind != -1)
					this._cursorPosition = ind+1;
			}
		});
		document.addEventListener("keydown", (e)=>{
			if (!this._isInTarget) return;
			// console.log(e);
			if (e.code!=e.key && !e.altKey && !e.ctrlKey && e.key != "Shift") {
				this.InsertStringAt(e.key, this._cursorPosition);
				this._cursorPosition++;
				for (let i = 0; i < this.#onChangeEvents.length; i++) {
					this.#onChangeEvents[i](this.value);
				}
			}
			if (e.key == "Backspace" && this._cursorPosition > 0) {
				this.SetText(this.value.slice(0,this._cursorPosition-1)+this.value.slice(this._cursorPosition));
				this._cursorPosition--;
				for (let i = 0; i < this.#onChangeEvents.length; i++) {
					this.#onChangeEvents[i](this.value);
				}
			}
			if (e.key == "Delete" && this._cursorPosition < this.value.length) {
				this.SetText(this.value.slice(0,this._cursorPosition)+this.value.slice(this._cursorPosition+1));
				// this._cursorPosition--;
				for (let i = 0; i < this.#onChangeEvents.length; i++) {
					this.#onChangeEvents[i](this.value);
				}
			}
			if (e.key == "ArrowLeft") {
				this._cursorPosition--;
			}
			if (e.key == "ArrowRight") {
				this._cursorPosition++;
			}
			if (this._cursorPosition < 0) this._cursorPosition = 0;
			if (this._cursorPosition > this.value.length) this._cursorPosition = this.value.length;
		});
	}
	AddOnChangeEvent(func) {
		let ind = this.#onChangeEvents.push(func);
		return ind-1;
	}
	GetLocalPos(pos) {
		return [
			this.position[0]+pos[0],
			this.position[1]+pos[1]
		];
	}
	InsertStringAt(str,i) {
		this.SetText(this.value.slice(0,i)+str+this.value.slice(i));
	}
	SetText(value) {
		ctx.font = this.fontSize+"px Arial";
		this.value = value;
	}
	DrawBackground() {
		d.rect(this.position[0],this.position[1],this.width,this.height,"black","white",1,true);
		d.rect(this.position[0],this.position[1],this.width,this.height,"black","rgba(255,255,255,0.5)",2,false);
	}
	GetTextPixelWidth(txt = this.value) {
		ctx.font = this.fontSize+"px Arial";
		return ctx.measureText(txt).width;
	}
	DrawTextCursor() {
		let wid = this.GetTextPixelWidth(this.value.slice(0,this._cursorPosition));
		let linePos1 = this.GetLocalPos([wid,0]);
		let linePos2 = this.GetLocalPos([wid,this.height]);
		d.line(linePos1[0],linePos1[1]+2, linePos2[0],linePos2[1]-2, this._textCursorColor,1);
	}
	DrawText() {
		ctx.font = this.fontSize+"px Arial";
		ctx.fillStyle = "White";
		ctx.fillText(this.value,this.position[0]+1,this.position[1]+this.fontSize);
	}
	Draw() {
		this._pixelTextLength = this.GetTextPixelWidth(this.value);
		this.DrawBackground();
		if (this._isInTarget)
			this.DrawTextCursor();
		this.DrawText();
	}
}
class NewGraph {
	constructor(DrawPosition, Width, Height, ScaleX, ScaleY) {
		this.position = DrawPosition;
		this.size = [Width, Height];
		this.gridSize = 4;
		this.scale = [ScaleX, ScaleY];

		this.data = []; // [{color:"red", name:"Animals count", data: [], bounded: 100}, ...], data: [{x: 0, y: 1}, style: "line"]
		this.equations = []; // [{color:"Red",name:"sin wave", func: (x)=>{}, style: "line"}];
		this.vectorField = undefined; // (x,y) => {return [x,y]; };
		this.colorField = undefined; // (x,y) => {return [r,g,b]; };

		this.centere = [this.position[0] + this.size[0] / 2, this.position[1] + this.size[1] / 2];
		this.viewOffset = [0,0];

		document.addEventListener("wheel", e => {
			this.OnScroll(e);
		});
		
		document.addEventListener("mousemove", e => {
			this.OnMouseDrag(e);
		});
	}
	SetPosition(p) {
		this.position = p;
		this.centere = [this.position[0] + this.size[0] / 2, this.position[1] + this.size[1] / 2];
	}


	// View offset exist in canvas space
	TransformGraphToFunc(p) {
		return [
			(p[0]-this.size[0]/2)*2*this.scale[0]/this.size[0],
			(p[1]-this.size[1]/2)*2*this.scale[1]/this.size[1],
		];
	}
	TransformFuncToGraph(p) {
		return [
			this.centere[0]+p[0]*this.size[0]/2/this.scale[0] - this.position[0],
			this.centere[1]-p[1]*this.size[1]/2/this.scale[1] - this.position[1],
		];
	}
	TransformCanvasToGraph(p) {
		return [
			p[0] - this.position[0],
			p[1] - this.position[1],
		];
	}
	TransformGraphToCanvas(p, sumView = false) {
		if (sumView) {
			p[0] += this.viewOffset[0];
			p[1] += this.viewOffset[1];
		}
		return [
			p[0]+this.position[0],
			p[1]+this.position[1],
		];
	}
	TransformToGlobalCanvas(p) {
		if (typeof (p) == "number") {
			return (p*this.size[0]/this.scale[0]/2)// / (this.size[0]/this.size[1]) * (this.scale[0]/this.scale[1]);
		}
		return [
			this.centere[0] + p[0]*this.size[0]/this.scale[0]/2 + this.viewOffset[0],
			this.centere[1] - p[1]*this.size[1]/this.scale[1]/2 + this.viewOffset[1],
		];
	}
	TransformToGlobalCanvasDir(p) {
		return [
			p[0]*this.size[0]/this.scale[0]/2,
			-p[1]*this.size[1]/this.scale[1]/2,
		]; 
	}
	TransformToGraphRef(p) {
		return [
			(p[0] - this.centere[0] - this.viewOffset[0]) / (graph.size[0]/graph.scale[0]/2),
			(p[1] - this.centere[1] - this.viewOffset[1]) / (graph.size[1]/graph.scale[1]/2)
		]
	}


	SetViewOffset(x,y) {
		this.viewOffset = [
			-this.size[0]/this.scale[0]*x/2,
			-this.size[1]/this.scale[1]*y/2
		];
	}
	GetViewOffset() {
		return [
			-this.scale[0]/this.size[0]*this.viewOffset[0]*2,
			-this.scale[1]/this.size[1]*this.viewOffset[1]*2
		];
	}
	DrawBackground() {
		d.rect(this.position[0], this.position[1], this.size[0], this.size[1], "black");
		d.rect(this.position[0], this.position[1], this.size[0], this.size[1], "black", "rgba(255,255,255,0.4)", 1, false);
	}
	DrawHUD() {
		for (let i = 0; i < this.data.length; i++) {

		}
	}
	DrawAxisLines() {
		// Draw in functional space
		let startPosition = this.TransformGraphToCanvas(this.TransformFuncToGraph([0,0]));

		let positiveBorderPosition = this.TransformGraphToCanvas(this.TransformFuncToGraph([this.scale[0],this.scale[1]]));
		let negativeBorderPosition = this.TransformGraphToCanvas(this.TransformFuncToGraph([-this.scale[0],-this.scale[1]]));
		d.ray(negativeBorderPosition[0], startPosition[1]+this.viewOffset[1], { x: 1, y: 0 }, positiveBorderPosition[0]-negativeBorderPosition[0], "rgba(255,255,255,0.4)", true, 0.95, 0.015);
		d.ray(startPosition[0]+this.viewOffset[0], negativeBorderPosition[1], { x: 0, y: 1}, positiveBorderPosition[1]-negativeBorderPosition[1], "rgba(255,255,255,0.4)", true, 0.95, 0.015);

		d.txt(this.scale[0],positiveBorderPosition[0],startPosition[1],"16px Arial","white");
		d.txt(this.scale[1],startPosition[0],positiveBorderPosition[1],"16px Arial","white");

		let viewRepeat = [this.size[0]/this.gridSize/2,this.size[1]/this.gridSize/2];
		ctx.font = "12px Arial";

		for (let x0 = -1; x0 <= 1; x0 += 1/this.gridSize) {
			let x = x0*this.scale[0];
			let xl = this.TransformGraphToCanvas(this.TransformFuncToGraph([x,0]))[0];
			if (xl+this.viewOffset[0]%viewRepeat[0] < this.position[0] || xl+this.viewOffset[0]%viewRepeat[0] > this.position[0]+this.size[0]) continue;

			// lines
			d.line(xl+this.viewOffset[0]%viewRepeat[0], this.position[1], xl+this.viewOffset[0]%viewRepeat[0], this.position[1]+this.size[1], "rgba(255,255,255,0.15)");

			// txt value
			let txtY = startPosition[1]+this.viewOffset[1];
			if (txtY < this.position[1]) txtY = this.position[1];
			if (txtY > this.position[1]+this.size[1]) txtY = this.position[1]+this.size[1];
			
			let txtValueOffset = -Math.floor(this.viewOffset[0]/this.size[0]*this.gridSize*2)/this.gridSize*this.scale[0];
			if (this.viewOffset[0] < 0) txtValueOffset -= 1/this.gridSize*this.scale[0];

			if (x+txtValueOffset != 0)
				d.txt((x+txtValueOffset).toFixed(2), xl+this.viewOffset[0]%viewRepeat[0]-ctx.measureText((x+txtValueOffset).toFixed(2)).width/2, txtY,"","rgba(255,255,255,0.6)");
		}
		for (let y0 = -1; y0 <= 1; y0 += 1/this.gridSize) {
			let y = y0*this.scale[1];
			let yl = this.TransformGraphToCanvas(this.TransformFuncToGraph([0,y]))[1];
			if (yl+this.viewOffset[1]%viewRepeat[1] <= this.position[1] || yl+this.viewOffset[1]%viewRepeat[1] >= this.position[1]+this.size[1]) continue;

			// lines
			d.line(this.position[0], yl+this.viewOffset[1]%viewRepeat[1], this.position[0]+this.size[0], yl+this.viewOffset[1]%viewRepeat[1], "rgba(255,255,255,0.15)");
			
			// txt value
			let txtX = startPosition[0]+this.viewOffset[0];
			if (txtX < this.position[0]) txtX = this.position[0];
			if (txtX > this.position[0]+this.size[0]) txtX = this.position[0]+this.size[0];

			let txtValueOffset = Math.floor(this.viewOffset[1]/this.size[1]*this.gridSize*2)/this.gridSize*this.scale[1];
			if (this.viewOffset[1] < 0) txtValueOffset += 1/this.gridSize*this.scale[1];

			if (y+txtValueOffset != 0)
				d.txt((y+txtValueOffset).toFixed(2), txtX, yl+this.viewOffset[1]%viewRepeat[1],"","rgba(255,255,255,0.6)");
		}
	}
	DrawDataAboutPoint(data,x,y, name = "") {
		let clr = "rgba(255,255,255,0.75)";

		d.txt(""+name,x,y,"12px Arial",clr);

		let pos = this.TransformGraphToCanvas(this.TransformFuncToGraph([data.x,data.y]),true);
		let leftPos = this.TransformGraphToCanvas(this.TransformFuncToGraph([0,data.y]),true);
		let downPos = this.TransformGraphToCanvas(this.TransformFuncToGraph([data.x,0]),true);
		if (leftPos[0] < this.position[0]) leftPos[0] = this.position[0];
		if (downPos[1] > this.position[1]+this.size[1]) downPos[1] = this.position[1]+this.size[1];

		ctx.strokeStyle = clr;
		ctx.beginPath();
		ctx.setLineDash([5, 5]);
		
		ctx.moveTo(leftPos[0],leftPos[1]);
		ctx.lineTo(pos[0],pos[1]);
		ctx.lineTo(downPos[0],downPos[1]);

		ctx.stroke();
		ctx.setLineDash([0,0]);
		

		let xTxtPos = (leftPos[0]+pos[0])/2 - ctx.measureText("X: " + data.x.toFixed(3)).width/2;
		let yTxtPos = (downPos[1]+pos[1])/2;
		d.txt("X: "+data.x.toFixed(3), xTxtPos, leftPos[1],"12px Arial",clr);
		d.txt("Y: "+data.y.toFixed(3), downPos[0], yTxtPos,"12px Arial",clr);
		// d.txt("Y: "+data.y,,,"12px Arial","white");
	}
	DrawDataAboutEquation(data, maxInd) {
		let area = 0;
		for (let i = 0; i < maxInd; i++)  {
			area += data.data[i].y;
		}
		area /= maxInd;
		// d.txt("Area size: " + area.toFixed(-Math.floor(Math.log10(area))+1).toString(),this.position[0]+this.size[0],this.position[1],"16px Arial", "white");
	}
	DrawData() {
		let dataOfMin= {};
		let nameOfMin = "";
		let posOfMin = [0,0];
		for (let D in this.data) {
			let DD = this.data[D];
			if (typeof(DD.bounded) != 'undefined' && DD.data.length > DD.bounded-1) {
				const sDelta = (DD.data[1].x-DD.data[0].x);
				const lenDelta = Math.abs(DD.bounded-DD.data.length);
				DD.data.splice(0,lenDelta);
				// console.log((DD.data[1].x-DD.data[0].x));
				for (let i = 0; i < DD.data.length; i++) {
					DD.data[i].x -= sDelta*(lenDelta);
				}
			}

			let val = this.data[D];
			let minDst = Infinity;
			let dataInd = -1;
			for (let j = 1; j < val.data.length; j++) {
				let p1 = this.TransformGraphToCanvas(this.TransformFuncToGraph([val.data[j-1].x,val.data[j-1].y]),true);
				let p2 = this.TransformGraphToCanvas(this.TransformFuncToGraph([val.data[j].x,val.data[j].y]),true);

				let dst = math.distanceToLine([input.mouse.x,input.mouse.y], p1, p2);
				dst[1] /= math.magnitude([p2[0]-p1[0],p2[1]-p1[1]]);
				if (dst[0] < 5 && minDst > dst[0]) {
					minDst = dst[0];
					dataOfMin =	{
						x: math.lerp(val.data[j-1].x,val.data[j].x,dst[1]),
						y: math.lerp(val.data[j-1].y,val.data[j].y,dst[1])
					};
					nameOfMin = val.name;
					dataInd = j-1;
					// posOfMin[0]= (p1[0]+p2[0])/2;
					// posOfMin[1]= (p1[1]+p2[1])/2;
					posOfMin[0]= input.mouse.x+12;
					posOfMin[1]= input.mouse.y;
				}
			}
			if (val.style == "dot") {
				for (let j = 0; j < val.data.length; j++) {
					let p1 = this.TransformGraphToCanvas(this.TransformFuncToGraph([val.data[j].x,val.data[j].y]),true);
	
					if ((Math.abs(p1[0]-this.centere[0]) > this.size[0]/2) ||
					(Math.abs(p1[1]-this.centere[1]) > this.size[1]/2)) continue;
				
					ctx.strokeStyle = "yellow";
					if (j > dataInd) ctx.strokeStyle = val.color;
					
					ctx.strokeRect(p1[0]-1.5,p1[1]-1.5,3,3);
					ctx.beginPath();
					ctx.moveTo(p1[0],this.centere[1]+this.viewOffset[1]);
					ctx.lineTo(p1[0],p1[1]);
					ctx.stroke();
				}
			} else {
				for (let j = 1; j < val.data.length; j++) {
					let p1 = this.TransformGraphToCanvas(this.TransformFuncToGraph([val.data[j-1].x,val.data[j-1].y]),true);
					let p2 = this.TransformGraphToCanvas(this.TransformFuncToGraph([val.data[j].x,val.data[j].y]),true);
	
					if ((Math.abs(p1[0]-this.centere[0]) > this.size[0]/2 && Math.abs(p2[0]-this.centere[0]) > this.size[0]/2) ||
					(Math.abs(p1[1]-this.centere[1]) > this.size[1]/2 && Math.abs(p2[1]-this.centere[1]) > this.size[1]/2)) continue;
				
					ctx.strokeStyle = "yellow";
					if (j > dataInd) ctx.strokeStyle = val.color;
					
					ctx.beginPath();
					ctx.moveTo(p1[0],p1[1]);
					ctx.lineTo(p2[0],p2[1]);
					ctx.stroke();
				}
			}
			/*for (let j = 0; j < this.data[i].data.length; j++) {
				let p = this.LocalToGlobal(this.data[i].data[j].x,this.data[i].data[j].y);
				// p[0] += this.viewOffset[0];
				// p[1] += this.viewOffset[1];
				if ((Math.abs(p[0]-this.centere[0]) > this.size[0]/2) || (Math.abs(p[1]-this.centere[1]) > this.size[1]/2)) continue;
				if (dataInd == j || dataInd-1 == j) 
					ctx.strokeStyle = "Blue";
				else
					ctx.strokeStyle = this.data[i].color;
				ctx.beginPath();
				ctx.arc(p[0],p[1],0.1/this.scale[1],0,2*Math.PI);
				ctx.stroke();
			}*/
			if (minDst < Infinity) {
				this.DrawDataAboutPoint(dataOfMin, posOfMin[0],posOfMin[1],nameOfMin);
				this.DrawDataAboutEquation(val, dataInd);
			}
		}
	}
	DrawColorField() {
		// if (this.scale[0] > 60 || this.scale[1] > 60) return;
		// let blockSize= [
		// 	this.size[0]/this.scale[0],
		// 	this.size[1]/this.scale[1]
		// ];
		// for (let y = -1; y <= 1; y+= 1/this.scale[1]) {
		// 	for (let x = -1; x <= 1; x+= 1/this.scale[0]) {
		// 		let v = this.colorField(x-1/this.scale[0]*0.5-this.viewOffset[0]/this.size[0]*2,-y+this.viewOffset[1]/this.size[1]*2);
		// 		if (v[0] < 0) v[0] = 0;
		// 		if (v[1] < 0) v[1] = 0;
		// 		if (v[2] < 0) v[2] = 0;
		// 		let drawPos = [
		// 			this.centere[0] + this.size[0]*x/2 + this.viewOffset[0]%(this.size[0]/this.scale[0]/2),
		// 			this.centere[1] + this.size[1]*y/2 + this.viewOffset[1]%(this.size[1]/this.scale[1]/2),
		// 		];
		// 		d.rect(drawPos[0]-blockSize[0]/2,drawPos[1]-blockSize[1]/2,blockSize[0],blockSize[1],inRgb(v[0]*255,v[1]*255,v[2]*255,255),"white");
		// 	}
		// }
		let fieldDensity = 30;
		for (let y = -1; y <= 1; y+=1/fieldDensity) {
			for (let x = -1; x <= 1; x+=1/fieldDensity) {
				let p = [x*this.scale[0]-this.viewOffset[0]/(this.size[0]/this.scale[0]/2),y*this.scale[1]-this.viewOffset[1]/(this.size[1]/this.scale[1]/2)];
				let v = this.colorField(p[0],p[1]);

				let drawPos = [this.centere[0]+x*(this.size[0]/2),this.centere[1]+y*(this.size[1]/2)];
				d.rect(drawPos[0],drawPos[1],this.size[0]/2/(fieldDensity)+1,this.size[1]/2/(fieldDensity)+1,inRgb(v[0]*255,v[1]*255,v[2]*255));
			}
		}
	}
	DrawVectorField() {
		// if (this.scale[0] > 60 || this.scale[1] > 60) return;
		let blockSize= [
			this.size[0]/this.scale[0],
			this.size[1]/this.scale[1]
		];
		// for (let y = -1; y <= 1; y+= 1/this.scale[1]) {
		// 	for (let x = -1; x <= 1; x+= 1/this.scale[0]) {
		// 		let v = this.vectorField(x-this.viewOffset[0]/this.size[0]*2,-y+this.viewOffset[1]/this.size[1]*2);
		// 		let drawPos = [
		// 			this.centere[0] + this.size[0]*x/2 + this.viewOffset[0]%(this.size[0]/this.scale[0]/2),
		// 			this.centere[1] + this.size[1]*y/2 + this.viewOffset[1]%(this.size[1]/this.scale[1]/2),
		// 		];
		// 		d.ray(drawPos[0],drawPos[1],{x:v[0],y:v[1]},math.magnitude(v),"red");
		// 	}
		// }
		
		let fieldDensity = 10;
		for (let y = -1; y <= 1; y+= 1/fieldDensity) {
			for (let x = -1; x <= 1; x+= 1/fieldDensity) {
				let v = this.vectorField(x-this.viewOffset[0]/this.size[0]*2,-y+this.viewOffset[1]/this.size[1]*2);
				let drawPos = [
					this.centere[0] + this.size[0]*x/2 + this.viewOffset[0]%(this.size[0]/fieldDensity/2),
					this.centere[1] + this.size[1]*y/2 + this.viewOffset[1]%(this.size[1]/fieldDensity/2),
				];
				d.ray(drawPos[0],drawPos[1],{x:v[0],y:-v[1]},math.magnitude(v),"red");
			}
		}
	}
	DrawEquations() {
		for (let i = 0; i < this.equations.length; i++) {
			let pastPos = [0,0];

			let minDst = Infinity;
			let minData = {};
			
			let pos = [-this.size[0]/2,this.equations[i].func(-this.scale[0]-this.viewOffset[0]/this.size[0]*this.scale[0]*2)*this.size[1]/this.scale[1]/2-this.viewOffset[1]];
			ctx.strokeStyle = this.equations[i].color;
			for (let x = -1; x <= 1; x+=0.01) {
				pastPos = [pos[0],pos[1]];
				let funcX = x*this.scale[0]-this.viewOffset[0]/this.size[0]*2*this.scale[0];
				pos = [x*this.size[0]/2+this.centere[0],-this.equations[i].func(funcX)*this.size[1]/this.scale[1]/2+this.viewOffset[1]+this.centere[1]];
				
				// if ((Math.abs(pos[0]-this.centere[0]) > this.size[0]/2) || (Math.abs(pos[1]-this.centere[1]) > this.size[1]/2)) continue;
				if (pos[1] <= this.position[1] || pos[1] >= this.position[1]+this.size[1] ||
					pos[0] <= this.position[0] || pos[0] >= this.position[0]+this.size[0]
				) {
					
					continue;
				}

				let dst = math.distanceToLine([input.mouse.x,input.mouse.y], pastPos, pos);

				if (dst[0] < 5 && dst[0] < minDst) {
					minDst = dst[0];
					minData = {x:x*this.scale[0]+this.GetViewOffset()[0], y:this.equations[i].func(x*this.scale[0]+this.GetViewOffset()[0])};
				}

				if (typeof (this.equations[i].color) != 'string') {
					let clr = this.equations[i].color(funcX);
					ctx.strokeStyle = `rgb(${clr[0]*255},${clr[1]*255},${clr[2]*255})`;
				}
				ctx.beginPath();
				ctx.moveTo(pastPos[0],pastPos[1]);
				ctx.lineTo(pos[0],pos[1]);
				ctx.stroke();
			}
			if (minDst < Infinity) {
				this.DrawDataAboutPoint(minData, input.mouse.x, input.mouse.y, this.equations[i].name);
			}
		}
	}
	Draw() {
		this.DrawBackground();
		this.DrawAxisLines();
		if (typeof(this.colorField) == "undefined" && typeof(this.vectorField) == "undefined") {
		} else if (typeof(this.colorField) != "undefined" && typeof(this.vectorField) == "undefined") {
			this.DrawColorField();
			this.DrawAxisLines();
		} else if (typeof(this.colorField) == "undefined" && typeof(this.vectorField) != "undefined") {
			this.DrawVectorField();
			this.DrawAxisLines();
		} else {
			console.error("Graph has vector field and color field in same time");
		}
		this.DrawData();
		this.DrawEquations();
	}
	OnScroll(e) {
		if (input.mouse.x < this.position[0] || input.mouse.x > this.position[0] + this.size[0]) return;
		if (input.mouse.y > this.position[1]+this.size[1] || input.mouse.y < this.position[1]) return;

		if (e.deltaY > 0) {
			this.viewOffset[0] /= 2;
			this.viewOffset[1] /= 2;

			this.scale[0] *= 2;
			this.scale[1] *= 2;
		} else if (e.deltaY < 0) {
			this.scale[0] *= 0.5;
			this.scale[1] *= 0.5;
			
			this.viewOffset[0] *= 2;
			this.viewOffset[1] *= 2;
		}
	}
	OnMouseDrag(e) {
		if (input.mouse.x < this.position[0] || input.mouse.x > this.position[0] + this.size[0]) return;
		if (input.mouse.y > this.position[1]+this.size[1] || input.mouse.y < this.position[1]) return;

		if (e.buttons == 4) {
			this.viewOffset[0] += e.movementX*0.5;
			this.viewOffset[1] += e.movementY*0.5;
		}
	}
	static WrapDataInBorders(Data,x1,x2, c = 100) {
		let newData = [];
		for (let i = 0; i < Data.length; i++) {
			newData.push({x:(x1+i/Data.length*x2),y:Data[i].y});
		}
		return newData;
	}
	static SliceData(Data,len,xStep = 1) {
		if (Data.length > len) {
			let c = Data.length-len;
			Data.splice(0,c);
			for (let i = 0; i < Data.length; i++) {
				Data[i].x-=xStep*c;
			}
		}
	}
	static ParseDataAsVector(data, xFunc = (x, len)=>{return x;}) { // data = [1,2,4,8]; return = [{x:1,y:1},{x:2,y:2},{x:3,y:4},{x:4,y:8}]
		let newData = [];
		for (let i = 0; i < data.length; i++) {
			newData.push({x:xFunc(i,data.length), y: data[i]});
		}
		return newData;
	}
	static ParseEquationAsData(func, dataLength, bounds = [0,1]) {
		let data = [];
		for (let i = 0; i < dataLength; i++) {
			let t = bounds[0]+i/(dataLength-1)*(bounds[1]-bounds[0]);
			data.push(func(t));
		}
		return data;
	}
	static ArrayToObject(arr) {
		let D = [];
		for (let i = 0; i < arr.length; i++) {
			D.push({x:arr[i][0],y:arr[i][1]});
		}
		return D;
	}
}
class Graph3d {
	constructor(cam) {
		this.position = [0, 0, 0];

		this.scale = [1, 1, 1];

		this.planeEquations = [];
		/* 
			{
				Name: name: "",
				Function: func: (x,y) => {},
				Color: clr: rgb(255,0,0)
			}
		*/

		this.gridPlane = [];
		/*
			{
				Name: name: "",
				Color: clr: rgb(255,0,0),
				Dots: grid: [[x1,x2,...],[x1,x2,...],...]
			}
		*/

		this.equationDrawQuality = 0.05;
		this.camera = cam;

		document.addEventListener("wheel", (e) => {
			let v = -e.deltaY / 100;
			this.scale[0] *= Math.pow(1.2, v);
			this.scale[1] *= Math.pow(1.2, v);
			this.scale[2] *= Math.pow(1.2, v);
		});
	}
	DrawAxes() {
		let xp1 = this.camera.ProjectToCanvas({ x: -1, y: 0, z: 0 });
		let xp2 = this.camera.ProjectToCanvas({ x: 1, y: 0, z: 0 });

		let yp1 = this.camera.ProjectToCanvas({ x: 0, y: -1, z: 0 });
		let yp2 = this.camera.ProjectToCanvas({ x: 0, y: 1, z: 0 });

		let zp1 = this.camera.ProjectToCanvas({ x: 0, y: 0, z: -1 });
		let zp2 = this.camera.ProjectToCanvas({ x: 0, y: 0, z: 1 });

		d.line(xp1.x, xp1.y, xp2.x, xp2.y, "red");
		d.line(yp1.x, yp1.y, yp2.x, yp2.y, "green");
		d.line(zp1.x, zp1.y, zp2.x, zp2.y, "blue");

		d.txt("X " + this.scale[0].toFixed(2), xp2.x, xp2.y, (30 / this.camera.DistToCamera({ x: 1, y: 0, z: 0 })) + "px Arial", "red");
		d.txt("Y " + this.scale[1].toFixed(2), yp2.x, yp2.y, (30 / this.camera.DistToCamera({ x: 0, y: 1, z: 0 })) + "px Arial", "green");
		d.txt("Z " + this.scale[2].toFixed(2), zp2.x, zp2.y, (38 / this.camera.DistToCamera({ x: 0, y: 0, z: 1 })) + "px Arial", "blue");
	}
	DrawPlaneEquation() {
		for (let i = 0; i < this.planeEquations.length; i++) {
			ctx.strokeStyle = this.planeEquations[i].clr;
			ctx.beginPath();
			for (let y = -this.scale[1]; y < this.scale[1]; y += this.equationDrawQuality * this.scale[1]) {
				let vp0 = this.camera.ProjectToCanvas({ x: -1, y: this.planeEquations[i].func(-this.scale[0], y) / this.scale[1], z: y / this.scale[2] });
				ctx.moveTo(vp0.x, vp0.y);
				for (let x = -this.scale[0]; x < this.scale[0]; x += this.equationDrawQuality * this.scale[0]) {
					let z = this.planeEquations[i].func(x, y);

					let vp = this.camera.ProjectToCanvas({ x: x / this.scale[0], y: z / this.scale[1], z: y / this.scale[2] });
					ctx.lineTo(vp.x, vp.y);
				}
			}
			for (let y = -this.scale[1]; y < this.scale[1]; y += this.equationDrawQuality * this.scale[1]) {
				let vp0 = this.camera.ProjectToCanvas({ x: y / this.scale[2], y: this.planeEquations[i].func(y, -this.scale[0]) / this.scale[1], z: -1 });
				ctx.moveTo(vp0.x, vp0.y);
				for (let x = -this.scale[0]; x < this.scale[0]; x += this.equationDrawQuality * this.scale[0]) {
					let z = this.planeEquations[i].func(y, x);

					let vp = this.camera.ProjectToCanvas({ x: y / this.scale[0], y: z / this.scale[1], z: x / this.scale[2] });
					ctx.lineTo(vp.x, vp.y);
				}
			}
			ctx.stroke();
		}

	}
	DrawPlaneGrid() {
		for (let i = 0; i < this.gridPlane.length; i++) {
			ctx.strokeStyle = this.gridPlane[i].clr;
			ctx.beginPath();
			for (let j = 0; j < this.gridPlane[i].grid.length; j++) {
				let vp0 = this.camera.ProjectToCanvas({ x: (j / this.gridPlane[i].grid.length - 0.5) * 2, y: (this.gridPlane[i].grid[j][0]), z: -1 });
				ctx.moveTo(vp0.x, vp0.y);
				for (let k = 0; k < this.gridPlane[i].grid[j].length; k++) {

					let vp = this.camera.ProjectToCanvas({ x: (j / this.gridPlane[i].grid.length - 0.5) * 2, y: this.gridPlane[i].grid[j][k], z: (k / this.gridPlane[i].grid[j].length - 0.5) * 2 });
					ctx.lineTo(vp.x, vp.y);
				}
			}
			for (let j = 0; j < this.gridPlane[i].grid.length; j++) {
				let vp0 = this.camera.ProjectToCanvas({ x: -1, y: (this.gridPlane[i].grid[0][j]), z: (j / this.gridPlane[i].grid.length - 0.5) * 2 });
				ctx.moveTo(vp0.x, vp0.y);
				for (let k = 0; k < this.gridPlane[i].grid[j].length; k++) {

					let vp = this.camera.ProjectToCanvas({ x: (k / this.gridPlane[i].grid[j].length - 0.5) * 2, y: this.gridPlane[i].grid[k][j], z: (j / this.gridPlane[i].grid.length - 0.5) * 2 });
					ctx.lineTo(vp.x, vp.y);
				}
			}
			ctx.stroke();
		}

	}
}

function MatrixMult3x3(v, m) {
	return {
		x: v.x * m[0][0] + v.y * m[1][0] + v.z * m[2][0],
		y: v.x * m[0][1] + v.y * m[1][1] + v.z * m[2][1],
		z: v.x * m[0][2] + v.y * m[1][2] + v.z * m[2][2]
	};
}
class Camera {
	constructor(StartOffset = [0,0], StartZoom = 1, MoveButton = 1) {
		this.offset = StartOffset;
		this.zoom = StartZoom;

		this.isMoveDependsOnZoom = true;
		this.isControlLocked = false;

		this.MoveButton = MoveButton;

		document.addEventListener("mousemove", (e) => { this.OnMouseMove(e) });
		document.addEventListener("wheel", (e) => { this.OnMouseScroll(e) });
	}
	GlobalToScreen(p) {
		if (typeof p.x == "number") {
			return {
				x: ((p.x - canv.width / 2 + this.offset[0]) * this.zoom) + canv.width / 2,
				y: ((p.y - canv.height / 2 + this.offset[1]) * this.zoom) + canv.height / 2
			};
		}
		if (Array.isArray(p)) {
			return [
				((p[0] - canv.width / 2 + this.offset[0]) * this.zoom) + canv.width / 2,
				((p[1] - canv.height / 2 + this.offset[1]) * this.zoom) + canv.height / 2
			];
		}
		if (typeof p == "number") {
			return p * this.zoom;
		}
		return null;
	}
	ScreenToGlobal(p) {
		if (typeof p.x == "number") {
			return {
				x: ((p.x - canv.width / 2 - this.offset[0]) / this.zoom - this.offset[0]) + canv.width / 2,
				y: ((p.y - canv.height / 2 - this.offset[1]) / this.zoom - this.offset[1]) + canv.height / 2
			};
		}
		if (Array.isArray(p)) {
			return [
				((p[0] - canv.width / 2 - this.offset[0]) / this.zoom - this.offset[0]) + canv.width / 2,
				((p[1] - canv.height / 2 - this.offset[1]) / this.zoom - this.offset[1]) + canv.height / 2
			];
		}
		if (typeof p == "number") {
			return p / this.zoom;
		}
	}
	OnMouseMove(e) {
		if (input.keyboard.code == 16) this.isControlLocked = true;
		else this.isControlLocked = false;
		if (e.buttons == this.MoveButton && !this.isControlLocked) {
			this.offset[0] += e.movementX / (this.isMoveDependsOnZoom ? this.zoom : 1);
			this.offset[1] += e.movementY / (this.isMoveDependsOnZoom ? this.zoom : 1);
		}
	}
	OnMouseScroll(e) {
		// if (input.keyboard.code == 16) this.isControlLocked = true;
		// else this.isControlLocked = false;
		// if (!this.isControlLocked) return;
		this.zoom += -e.deltaY / 100 * this.zoom * 0.1;
		if (this.zoom < 0.0001) this.zoom = 0.0001;
	}
}
class Camera3d {
	constructor(StartPosition = [0,0,0], StartFov = 0.8) {
		this.position = StartPosition;
		this.rotation = [0,0,0];
		this.fov = StartFov;

		this.forward = [0,0,1];
		this.right = [1,0,0];
		this.up = [0,1,0];

		document.addEventListener("mousemove", (e) => {
			this.Control(e.movementX,e.movementY);
		});
		// document.addEventListener();
	}
	Control(x,y) {
		if (input.mouse.click == 3) {
			this.rotation[0] += x * 0.003;
			this.rotation[1] += -y * 0.003;

			this.Rotate(this.rotation[0], this.rotation[1], this.rotation[2]);
		}
	}
	KeyboardControl(speed = 1) {
		if (input.keyboard.char == 'w') {
			this.position[0] += this.forward[0] * speed;
			this.position[1] += this.forward[1] * speed;
			this.position[2] += this.forward[2] * speed;
		}
		if (input.keyboard.char == 'a') {
			this.position[0] += -this.right[0] * speed;
			this.position[1] += -this.right[1] * speed;
			this.position[2] += -this.right[2] * speed;
		}
		if (input.keyboard.char == 's') {
			this.position[0] += -this.forward[0] * speed;
			this.position[1] += -this.forward[1] * speed;
			this.position[2] += -this.forward[2] * speed;
		}
		if (input.keyboard.char == 'd') {
			this.position[0] += this.right[0] * speed;
			this.position[1] += this.right[1] * speed;
			this.position[2] += this.right[2] * speed;
		}
		if (input.keyboard.char == 'r') {
			this.position[0] += this.up[0] * speed;
			this.position[1] += this.up[1] * speed;
			this.position[2] += this.up[2] * speed;
		}
		if (input.keyboard.char == 'f') {
			this.position[0] += -this.up[0] * speed;
			this.position[1] += -this.up[1] * speed;
			this.position[2] += -this.up[2] * speed;
		}
	}
	Rotate(x, y, z) {
		this.rotation[0] = x;
		this.rotation[1] = y;
		this.rotation[2] = z;

		// this -> forward = Vector3ff(cosf(rotation[1]) * sinf(rotation[0]), sinf(rotation[1]), cosf(rotation[0]) * cosf(rotation[1]));
		// this -> right = Vector3ff(sinf(rotation[0] + PI / 2.0f), 0.0f, cosf(rotation[0] + PI / 2.0f));
		// this -> up = Vector3ff(cosf(rotation[1] + PI / 2.0f) * sinf(rotation[0]), sinf(rotation[1] + PI / 2.0f), cosf(rotation[0]) * cosf(rotation[1] + PI / 2.0f));

		this.forward = [
			Math.cos(this.rotation[1]) * Math.sin(this.rotation[0]),
			Math.sin(this.rotation[1]),
			Math.cos(this.rotation[0]) * Math.cos(this.rotation[1])
		];
		this.right = [
			Math.sin(this.rotation[0] + Math.PI / 2),
			0,
			Math.cos(this.rotation[0] + Math.PI / 2)
		];
		this.up = [
			Math.cos(this.rotation[1] + Math.PI / 2) * Math.sin(this.rotation[0]),
			Math.sin(this.rotation[1] + Math.PI / 2),
			Math.cos(this.rotation[0]) * Math.cos(this.rotation[1] + Math.PI / 2)
		];
		return this.forward;

		/*let forward = {
			x: 0,
			y: 0,
			z: 1
		};
		this.forward = MatrixMult3x3(forward, [
			[1, 0, 0],
			[0, Math.cos(x), -Math.sin(x)],
			[0, Math.sin(x), Math.cos(x)],
		]);
		this.forward = MatrixMult3x3(this.forward, [
			[Math.cos(y), 0, Math.sin(y)],
			[0, 1, 0],
			[-Math.sin(y), 0, Math.cos(y)],
		]);
		this.forward = MatrixMult3x3(this.forward, [
			[Math.cos(z), -Math.sin(z), 0],
			[Math.sin(z), Math.cos(z), 0],
			[0, 0, 1],
		]);

		let right = {
			x: 1,
			y: 0,
			z: 0
		};
		this.right = MatrixMult3x3(right, [
			[1, 0, 0],
			[0, Math.cos(x), -Math.sin(x)],
			[0, Math.sin(x), Math.cos(x)],
		]);
		this.right = MatrixMult3x3(this.right, [
			[Math.cos(y), 0, Math.sin(y)],
			[0, 1, 0],
			[-Math.sin(y), 0, Math.cos(y)],
		]);
		this.right = MatrixMult3x3(this.right, [
			[Math.cos(z), -Math.sin(z), 0],
			[Math.sin(z), Math.cos(z), 0],
			[0, 0, 1],
		]);

		let up = {
			x: 0,
			y: 1,
			z: 0
		};
		this.up = MatrixMult3x3(up, [
			[1, 0, 0],
			[0, Math.cos(x), -Math.sin(x)],
			[0, Math.sin(x), Math.cos(x)],
		]);
		this.up = MatrixMult3x3(this.up, [
			[Math.cos(y), 0, Math.sin(y)],
			[0, 1, 0],
			[-Math.sin(y), 0, Math.cos(y)],
		]);
		this.up = MatrixMult3x3(this.up, [
			[Math.cos(z), -Math.sin(z), 0],
			[Math.sin(z), Math.cos(z), 0],
			[0, 0, 1],
		]);

		return this.forward;*/
	}
	isOnScreen(p) {
		if (p[0] < 0 || p[1] < 0) return false;
		if (p[0] > canv.width || p[1] > canv.height) return false;
		return true;
	}
	ProjectToCanvas(p) {
		let size = Math.min(canv.width, canv.height);

		let delta = [
			p[0] - this.position[0],
			p[1] - this.position[1],
			p[2] - this.position[2]
		];
		let nDelta = math.normalize(delta);
		if (math.dot(nDelta, this.forward) < 0) return [-1,-1];
		let v = [
			(nDelta[0] / math.dot(nDelta, this.forward) - this.forward[0]) * size * this.fov,
			(nDelta[1] / math.dot(nDelta, this.forward) - this.forward[1]) * size * this.fov,
			(nDelta[2] / math.dot(nDelta, this.forward) - this.forward[2]) * size * this.fov
		];

		return [
			math.dot(v, this.right) + canv.width / 2,
			-math.dot(v, this.up) + canv.height / 2
		];
	}
	OrtogonalProjectToCanvas(p) {
		let size = Math.min(canv.width, canv.height);
		let delta = [
			p[0] - this.position[0],
			p[1] - this.position[1],
			p[2] - this.position[2]
		];
		if (math.dot(delta, this.forward) < 0) return [-1,-1];
		return [
			math.dot(delta, this.right) * this.fov*size + canv.width / 2,
			-math.dot(delta, this.up) * this.fov*size + canv.height / 2
		];
	}
	RadialProjectToCanvas(p) {
		let size = Math.min(canv.width, canv.height);
		let delta = [
			p[0] - this.position[0],
			p[1] - this.position[1],
			p[2] - this.position[2]
		];
		let projP = math.normalize(delta);
		
		let deltaProj = [
			this.position[0]+this.forward[0]-projP[0],
			this.position[1]+this.forward[1]-projP[1],
			this.position[2]+this.forward[2]-projP[2]
		];

		return [-math.dot(deltaProj,this.right)*this.fov*size+canv.width/2,math.dot(deltaProj,this.up)*this.fov*size+canv.height/2];
	}
	ParabolicProjToCanvas(p, Attractor = [0,0]) {
		let size = Math.min(canv.width, canv.height);

		let delta = [
			p[0] - this.position[0],
			p[1] - this.position[1],
			p[2] - this.position[2]
		];
		let nDelta = math.normalize(delta);
		if (math.dot(nDelta, this.forward) < 0) return [-1,-1];
		
		let attProj = math.dot(this.forward,Attractor);
		let dirProj = math.dot(this.forward,nDelta);

		let D = attProj**2 + 4*dirProj;
		let t = -(attProj-Math.sqrt(D))/(2*dirProj);

		let P = [
			(nDelta[0]*t + Attractor[0])*t * size * this.fov,
			(nDelta[1]*t + Attractor[1])*t * size * this.fov,
		];

		return [
			math.dot(P,this.right) + canv.width/2,
			-math.dot(P,this.up) + canv.height/2
		];
	}
	ViewportPosToRay(p) {
		let size = Math.min(canv.width, canv.height);
		let deltas = {
			x: (p[0] - canv.width / 2) / size,
			y: -(p[1] - canv.height / 2) / size
		};
		let direction = {
			x: this.forward[0] + this.right[0] * deltas[0] + this.up[0] * deltas[1],
			y: this.forward[1] + this.right[1] * deltas[0] + this.up[1] * deltas[1],
			z: this.forward[2] + this.right[2] * deltas[0] + this.up[2] * deltas[1]
		};
		return direction;
	}
	DistToCamera(p) {
		return math.magnitude({
			x: p[0] - this.position[0],
			y: p[1] - this.position[1],
			z: p[2] - this.position[2]
		});
	}
	AxisLines() {
		let vp1 = this.ProjectToCanvas([-1,0,0]);
		let vp2 = this.ProjectToCanvas([1,0,0]);
		let vp3 = this.ProjectToCanvas([0,-1,0]);
		let vp4 = this.ProjectToCanvas([0,1,0]);
		let vp5 = this.ProjectToCanvas([0,0,-1]);
		let vp6 = this.ProjectToCanvas([0,0,1]);

		d.line(vp1[0], vp1[1], vp2[0], vp2[1], "red");
		d.line(vp3[0], vp3[1], vp4[0], vp4[1], "green");
		d.line(vp5[0], vp5[1], vp6[0], vp6[1], "blue");
	}
	DrawLine(p1,p2,clr, withCtx = false) {
		let proj1 = this.ProjectToCanvas(p1);
		let proj2 = this.ProjectToCanvas(p2);
		// if (this.isOnScreen(proj1) || this.isOnScreen(proj2))
		if (!withCtx)
			d.line(proj1[0],proj1[1], proj2[0],proj2[1], clr,1);
		else {
			ctx.moveTo(proj1[0],proj1[1]);
			ctx.lineTo(proj2[0],proj2[1]);
		}
	}
	DrawBox(p1,p2, clr = 'white', withCtx = false) {
		let corners = [
			[p1[0], p1[1], p1[2]],
			[p1[0], p1[1], p2[2]],
			[p1[0], p2[1], p2[2]],
			[p1[0], p2[1], p1[2]],
			[p1[0], p1[1], p1[2]],
			[p2[0], p1[1], p1[2]],
			[p2[0], p1[1], p2[2]],
			[p2[0], p2[1], p2[2]],
			[p2[0], p2[1], p1[2]],
			[p2[0], p1[1], p1[2]],
		];
		if (!withCtx) {
			let past = corners[0];
			for (let i = 1; i <= 4; i++) {
				cam.DrawLine(past,corners[i],clr, withCtx);
				past = corners[i];
			}
			past = corners[5];
			for (let i = 5; i <= 9; i++) {
				cam.DrawLine(past,corners[i],clr, withCtx);
				past = corners[i];
			}
			for (let i = 0; i < 4; i++) {
				cam.DrawLine(corners[i],corners[5+i],clr, withCtx);
			}
		}
		
	}
}
class Camera4d {
	constructor(StartPosition = { x: 0, y: 0, z: 0, w: 0 }, StartFov = 1) {
		this.position = StartPosition;
		this.rotation = { x: 0, y: 0, z: 0, w: 0 };
		this.fov = StartFov;

		this.forward = { x: 0, y: 0, z: 1, w: 0 };
		this.right = { x: 1, y: 0, z: 0, w: 0 };
		this.up = { x: 0, y: 1, z: 0, w: 0 };
	}
	Rotate(x, y, z, w) {
		this.rotation.x = x;
		this.rotation.y = y;
		this.rotation.z = z;
		let forward = {
			x: 0,
			y: 0,
			z: 1
		};
		this.forward = MatrixMult3x3(forward, [
			[1, 0, 0],
			[0, Math.cos(x), -Math.sin(x)],
			[0, Math.sin(x), Math.cos(x)],
		]);
		this.forward = MatrixMult3x3(this.forward, [
			[Math.cos(y), 0, Math.sin(y)],
			[0, 1, 0],
			[-Math.sin(y), 0, Math.cos(y)],
		]);
		this.forward = MatrixMult3x3(this.forward, [
			[Math.cos(z), -Math.sin(z), 0],
			[Math.sin(z), Math.cos(z), 0],
			[0, 0, 1],
		]);

		let right = {
			x: 1,
			y: 0,
			z: 0
		};
		this.right = MatrixMult3x3(right, [
			[1, 0, 0],
			[0, Math.cos(x), -Math.sin(x)],
			[0, Math.sin(x), Math.cos(x)],
		]);
		this.right = MatrixMult3x3(this.right, [
			[Math.cos(y), 0, Math.sin(y)],
			[0, 1, 0],
			[-Math.sin(y), 0, Math.cos(y)],
		]);
		this.right = MatrixMult3x3(this.right, [
			[Math.cos(z), -Math.sin(z), 0],
			[Math.sin(z), Math.cos(z), 0],
			[0, 0, 1],
		]);

		let up = {
			x: 0,
			y: 1,
			z: 0
		};
		this.up = MatrixMult3x3(up, [
			[1, 0, 0],
			[0, Math.cos(x), -Math.sin(x)],
			[0, Math.sin(x), Math.cos(x)],
		]);
		this.up = MatrixMult3x3(this.up, [
			[Math.cos(y), 0, Math.sin(y)],
			[0, 1, 0],
			[-Math.sin(y), 0, Math.cos(y)],
		]);
		this.up = MatrixMult3x3(this.up, [
			[Math.cos(z), -Math.sin(z), 0],
			[Math.sin(z), Math.cos(z), 0],
			[0, 0, 1],
		]);

		return this.forward;
	}
	ProjectToCanvas(p) {
		let size = Math.min(canv.width, canv.height);

		let delta = {
			x: p.x - this.position.x,
			y: p.y - this.position.y,
			z: p.z - this.position.z,
			w: p.w - this.position.w
		};
		let deltaLen = Math.sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z + delta.w * delta.w);
		let nDelta = {
			x: delta.x / deltaLen,
			y: delta.y / deltaLen,
			z: delta.z / deltaLen,
			w: delta.w / deltaLen
		};
		let t = (this.forward.x ** 2 + this.forward.y ** 2 + this.forward.z ** 2 + this.forward.w ** 2) / (nDelta.x * this.forward.x + nDelta.y * this.forward.y + nDelta.z * this.forward.z + nDelta.w * this.forward.w);
		let v = {
			x: (nDelta.x * t - this.forward.x) * size * this.fov,
			y: (nDelta.y * t - this.forward.y) * size * this.fov,
			z: (nDelta.z * t - this.forward.z) * size * this.fov,
			w: (nDelta.w * t - this.forward.w) * size * this.fov
		};
		return {
			x: (v.x * this.right.x + v.y * this.right.y + v.z * this.right.z + v.w * this.right.w) + canv.width / 2,
			y: -(v.x * this.up.x + v.y * this.up.y + v.z * this.up.z + v.w * this.up.w) + canv.height / 2
		};
		// return {
		// 	x: math.dot(v, this.right) + canv.width / 2,
		// 	y: -math.dot(v, this.up) + canv.height / 2
		// };
	}
	ViewportPosToRay(p) {
		let size = Math.min(canv.width, canv.height);
		let deltas = {
			x: (p.x - canv.width / 2) / size,
			y: -(p.y - canv.height / 2) / size
		};
		let direction = {
			x: this.forward.x + this.right.x * deltas.x + this.up.x * deltas.y,
			y: this.forward.y + this.right.y * deltas.x + this.up.y * deltas.y,
			z: this.forward.z + this.right.z * deltas.x + this.up.z * deltas.y
		};
		return direction;
	}
	DistToCamera(p) {
		return math.magnitude({
			x: p.x - this.position.x,
			y: p.y - this.position.y,
			z: p.z - this.position.z
		});
	}
}