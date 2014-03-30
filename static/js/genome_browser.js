var GenomeBrowser = function(htmlElementId, socket, width, height) {
    var self = this;
    this.canvas = Raphael(htmlElementId, width, height);
    this.genomeId = undefined;
    this.current_genome = undefined;
    this.socket = socket;
    this.genomicRange = [10000, 15000];
    this.genomicStep = (this.genomicRange[1] - this.genomicRange[0])/width;
    this.annotationsDisplayed = []; //a matrix (n,2) for the annotations displayed in the browser [[annotation, rect],[annotation, rect],...]
    this.annotation_to_blink = undefined;
    this.intervalBetweenTicks = Math.round((self.genomicRange[1]-self.genomicRange[0])/10);
    this.genomicSequence = undefined;
    this.expressionLevels = [];

    this.ticksBar = self.canvas.rect(0 , 15, width , 2);
    
    self.ticksBar.attr({
          fill: 'black',
          stoke: 0
        });
    
    var Tick = function(rank) {

        var label = self.canvas.text(rank/10*width, 35, "").attr({fill: 'black'});
      
        self.canvas.rect(-1+rank/10*width+(width/120)/2, 15, 2, 10).attr({
            fill:'black',
            stoke:0,
        });
     
        this.update = function() {
            var position = ""+Math.round(self.genomicRange[0]+rank/10*(self.genomicRange[1]-self.genomicRange[0]));
            if (position.length > 6) {
                if (self.intervalBetweenTicks >= 100000)
                    position = position.substring(0, position.length-6)+"."+position.substring(position.length-6, position.length-5)+"Mb";
                else if (self.intervalBetweenTicks >= 10000)
                    position = position.substring(0, position.length-6)+"."+position.substring(position.length-6, position.length-4)+"Mb";
                else if (self.intervalBetweenTicks >= 1000)
                    position = position.substring(0, position.length-6)+"."+position.substring(position.length-6, position.length-3)+"Mb";
                else if (self.intervalBetweenTicks >= 100)
                    position = position.substring(0, position.length-6)+"."+position.substring(position.length-6, position.length-2)+"Mb";
                else if (self.intervalBetweenTicks >= 12)
                    position = position.substring(0, position.length-6)+"."+position.substring(position.length-6, position.length)+"Mb";
            } else if (position.length > 3) {
                if (self.intervalBetweenTicks >= 100)
                    position = position.substring(0, position.length-3)+"."+position.substring(position.length-3, position.length-2)+"kb";
                else if (self.intervalBetweenTicks >= 12)
                    position = position.substring(0, position.length-3)+"."+position.substring(position.length-3, position.length)+"kb";
            }
            label.attr('text', position);   
        };

    };
    
    this.ticks = [new Tick(0), new Tick(1), new Tick(2), new Tick(3), new Tick(4), new Tick(5), new Tick(6), new Tick(7), new Tick(8), new Tick(9), new Tick(10)];

    this.ticks.forEach(function(tick) { //to get the first label
        tick.update();
    });

    this.residues = [];

    for (var i = 0 ; i < 120 ; i++) {
        var rect = self.canvas.rect(i*width/120 , 0, width/120 , 15);
        rect.attr({fill: 'red', opacity:0, 'stroke-width':0.1});
        var text = self.canvas.text(rect.getBBox().x+rect.getBBox().width/2, rect.getBBox().y2-rect.getBBox().height/2 , 'A').attr({fill: 'white', opacity:0});
        this.residues.push([text,rect]);    
    }

    var ExpressionLevel = function(rank) {
        
        var rect = self.canvas.rect(rank*10, height, 10, 0).attr({
            fill:'orange',
            "fill-opacity": .4, 
            stroke:'orange'
        });

        this.update = function(newValue) {
            rect.attr('y',height-50-newValue);
            rect.attr('height',newValue);
            rect.animate({"fill-opacity":0.4, "stroke-opacity":1}, 1000);
        };

        this.hide = function() {
            rect.attr({"fill-opacity":0, "stroke-opacity":0});
        };       
    };

    for (var i = 0 ; i < 120 ; i++)
        this.expressionLevels.push(new ExpressionLevel(i));

    this.ox = 0;
    this.oy = 0;

    this.header = 50; //where the ruler and the ticks will be drawn

    //the dock from top to down
    //first the border of the dock.
    this.canvas.rect(5,(height/2-130), 32 , 260).attr({stroke:'#000', fill:'white'});

    //go up
    this.canvas.rect(3, 0, 27, 30)
        .attr({fill: "#FFF", "stroke":"none"}).transform("t"+5+","+(height/2-125)).click(function() {
            self.goUpDown(-25);//25 = 20 px height for an annotation + 5px space between annotations
        });
    this.canvas.path('M25.682,24.316L15.5,6.684L5.318,24.316H25.682z')
        .attr({fill: "#000", stroke: "none", title : "Show top annotations"}).transform("t"+5+","+(height/2-125)).click(function() {
            self.goUpDown(-25);//25 = 20 px height for an annotation + 5px space between annotations
        });
    //large move left
    this.canvas.rect(3, 0, 27, 30)
        .attr({fill: "#FFF", "stroke":"none"}).transform("t"+5+","+(height/2-95)).click(function() {
            self.walkOnGenome(self.getCurrentGenomicCenter()-(self.intervalBetweenTicks*10));
        });
    this.canvas.path('M5.5,15.499,15.8,21.447,15.8,15.846,25.5,21.447,25.5,9.552,15.8,15.152,15.8,9.552z')
        .attr({fill: "#000", stroke: "none", title : "Faster move upstream"}).transform("t"+5+","+(height/2-95)).click(function() {
            self.walkOnGenome(self.getCurrentGenomicCenter()-(self.intervalBetweenTicks*10));
        });
    //small move left
    this.canvas.rect(3, 0, 27, 30)
        .attr({fill: "#FFF", "stroke":"none"}).transform("t"+5+","+(height/2-67)).click(function() {
            self.walkOnGenome(self.getCurrentGenomicCenter()-self.intervalBetweenTicks);
        });
    this.canvas.path('M24.316,5.318L6.684,15.5l17.632,10.182V5.318L24.316,5.318z')
        .attr({fill: "#000", stroke: "none", title : "Move upstream"}).transform("t"+5+","+(height/2-67)).click(function() {
            self.walkOnGenome(self.getCurrentGenomicCenter()-self.intervalBetweenTicks);
        });

    //zoom out
    this.canvas.rect(3, 0, 27, 30)
        .attr({fill: "#FFF", "stroke":"none"}).transform("t"+5+","+(height/2-35)).click(function() {
            self.scale(2);
        }); -
    this.canvas.path('M22.646,19.307c0.96-1.583,1.523-3.435,1.524-5.421C24.169,8.093,19.478,3.401,13.688,3.399C7.897,3.401,3.204,8.093,3.204,13.885c0,5.789,4.693,10.481,10.484,10.481c1.987,0,3.839-0.563,5.422-1.523l7.128,7.127l3.535-3.537L22.646,19.307zM13.688,20.369c-3.582-0.008-6.478-2.904-6.484-6.484c0.006-3.582,2.903-6.478,6.484-6.486c3.579,0.008,6.478,2.904,6.484,6.486C20.165,17.465,17.267,20.361,13.688,20.369zM8.854,11.884v4.001l9.665-0.001v-3.999L8.854,11.884z')
        .attr({fill: "#000", stroke: "none"}).transform("t"+5+","+(height/2-35)).click(function() {
            self.scale(2);
        });
    //zoom in
    this.canvas.rect(3, 0, 27, 30)
        .attr({fill: "#FFF", "stroke":"none"}).transform("t"+5+","+(height/2)).click(function() {
            self.scale(1/2);
        });
    this.canvas.path('M22.646,19.307c0.96-1.583,1.523-3.435,1.524-5.421C24.169,8.093,19.478,3.401,13.688,3.399C7.897,3.401,3.204,8.093,3.204,13.885c0,5.789,4.693,10.481,10.484,10.481c1.987,0,3.839-0.563,5.422-1.523l7.128,7.127l3.535-3.537L22.646,19.307zM13.688,20.369c-3.582-0.008-6.478-2.904-6.484-6.484c0.006-3.582,2.903-6.478,6.484-6.486c3.579,0.008,6.478,2.904,6.484,6.486C20.165,17.465,17.267,20.361,13.688,20.369zM15.687,9.051h-4v2.833H8.854v4.001h2.833v2.833h4v-2.834h2.832v-3.999h-2.833V9.051z')
        .attr({fill: "#000", stroke: "none"}).transform("t"+5+","+(height/2)).click(function() {
            self.scale(1/2);
        });
    //small move right
    this.canvas.rect(3, 0, 27, 30)
        .attr({fill: "#FFF", "stroke":"none"}).transform("t"+5+","+(height/2+35)).click(function() {
            self.walkOnGenome(self.getCurrentGenomicCenter()+self.intervalBetweenTicks);
        });
    this.canvas.path('M6.684,25.682L24.316,15.5L6.684,5.318V25.682z')
        .attr({fill: "#000", stroke: "none", title : "Move downstream"}).transform("t"+5+","+(height/2+35)).click(function() {
            self.walkOnGenome(self.getCurrentGenomicCenter()+self.intervalBetweenTicks);
        });
    //large move right
    this.canvas.rect(3, 0, 27, 30)
        .attr({fill: "#FFF", "stroke":"none"}).transform("t"+5+","+(height/2+65)).click(function() {
            self.walkOnGenome(self.getCurrentGenomicCenter()+(self.intervalBetweenTicks*10));
        });
    this.canvas.path('M25.5,15.5,15.2,9.552,15.2,15.153,5.5,9.552,5.5,21.447,15.2,15.847,15.2,21.447z')
        .attr({fill: "#000", stroke: "none", title : "Faster move downstream"}).transform("t"+5+","+(height/2+65)).click(function() {
            self.walkOnGenome(self.getCurrentGenomicCenter()+(self.intervalBetweenTicks*10));
        });
    //go down
    this.canvas.rect(3, 0, 27, 30)
        .attr({fill: "#FFF", "stroke":"none"}).transform("t"+5+","+(height/2+95)).click(function() {
            self.goUpDown(25);//25 = 20 px height for an annotation + 5px space between annotations
        });
    this.canvas.path('M5.318,6.684L15.5,24.316L25.682,6.684H5.318z')
        .attr({fill: "#000", stroke: "none", title : "Show bottom annotations"}).transform("t"+5+","+(height/2+95)).click(function() {
            self.goUpDown(25);//25 = 20 px height for an annotation + 5px space between annotations
        });

    //legend
    this.canvas.rect(25, height-25, 20, 20).attr({fill: '#0B6121', "fill-opacity": .4, stroke:'#0B6121'});
    this.canvas.text(65, height-15, "Gene").attr({fill: '#0B6121', 'font-size': 12, 'font-weight':'bold'});
    this.canvas.rect(85, height-25, 20, 20).attr({fill: '#088A08', "fill-opacity": .4, stroke:'#088A08'});
    this.canvas.text(125, height-15, "mRNA").attr({fill: '#088A08', 'font-size': 12, 'font-weight':'bold'});
    this.canvas.rect(150, height-25, 20, 20).attr({fill: '#A5DF00', "fill-opacity": .4, stroke:'#A5DF00'});
    this.canvas.text(185, height-15, "CDS").attr({fill: '#A5DF00', 'font-size': 12, 'font-weight':'bold'});
    this.canvas.rect(210, height-25, 20, 20).attr({fill: 'orange', "fill-opacity": .4, stroke:'orange'});
    this.canvas.text(250, height-15, "Intron").attr({fill: 'orange', 'font-size': 12, 'font-weight':'bold'});
    this.canvas.rect(275, height-25, 20, 20).attr({fill: 'purple', "fill-opacity": .4, stroke:'purple'});
    this.canvas.text(315, height-15, "tRNA").attr({fill: 'purple', 'font-size': 12, 'font-weight':'bold'});
    this.canvas.rect(335, height-25, 20, 20).attr({fill: 'red', "fill-opacity": .4, stroke:'red'});
    this.canvas.text(375, height-15, "rRNA").attr({fill: 'red', 'font-size': 12, 'font-weight':'bold'});
    this.canvas.rect(400, height-25, 20, 20).attr({fill: '#0040FF', "fill-opacity": .4, stroke:'#0040FF'});
    this.canvas.text(445, height-15, "ncRNA").attr({fill: '#0040FF', 'font-size': 12, 'font-weight':'bold'});
    this.canvas.rect(470, height-25, 20, 20).attr({fill: 'grey', "fill-opacity": .4, stroke:'grey'});
    this.canvas.text(510, height-15, "Misc").attr({fill: 'grey', 'font-size': 12, 'font-weight':'bold'});

    //powered by...
    this.canvas.text(width-60, height-15, "Powered by").attr({fill: "#000", 'font-size': 10, 'font-weight':'bold'});
    //the raphael logo
    this.canvas.path('M27.777,18.941c0.584-0.881,0.896-1.914,0.896-2.998c0-1.457-0.567-2.826-1.598-3.854l-6.91-6.911l-0.003,0.002c-0.985-0.988-2.35-1.6-3.851-1.6c-1.502,0-2.864,0.612-3.85,1.6H12.46l-6.911,6.911c-1.031,1.029-1.598,2.398-1.598,3.854c0,1.457,0.567,2.826,1.598,3.854l6.231,6.229c0.25,0.281,0.512,0.544,0.789,0.785c1.016,0.961,2.338,1.49,3.743,1.49c1.456,0,2.825-0.565,3.854-1.598l6.723-6.725c0.021-0.019,0.034-0.032,0.051-0.051l0.14-0.138c0.26-0.26,0.487-0.54,0.688-0.838c0.004-0.008,0.01-0.015,0.014-0.021L27.777,18.941zM26.658,15.946c0,0.678-0.197,1.326-0.561,1.879c-0.222,0.298-0.447,0.559-0.684,0.784L25.4,18.625c-1.105,1.052-2.354,1.35-3.414,1.35c-0.584,0-1.109-0.09-1.523-0.195c-2.422-0.608-5.056-2.692-6.261-5.732c0.649,0.274,1.362,0.426,2.11,0.426c2.811,0,5.129-2.141,5.415-4.877l3.924,3.925C26.301,14.167,26.658,15.029,26.658,15.946zM16.312,5.6c1.89,0,3.426,1.538,3.426,3.427c0,1.89-1.536,3.427-3.426,3.427c-1.889,0-3.426-1.537-3.426-3.427C12.886,7.138,14.423,5.6,16.312,5.6zM6.974,18.375c-0.649-0.648-1.007-1.512-1.007-2.429c0-0.917,0.357-1.78,1.007-2.428l2.655-2.656c-0.693,2.359-0.991,4.842-0.831,7.221c0.057,0.854,0.175,1.677,0.345,2.46L6.974,18.375zM11.514,11.592c0.583,4.562,4.195,9.066,8.455,10.143c0.693,0.179,1.375,0.265,2.033,0.265c0.01,0,0.02,0,0.027,0l-3.289,3.289c-0.648,0.646-1.512,1.006-2.428,1.006c-0.638,0-1.248-0.177-1.779-0.5l0.001-0.002c-0.209-0.142-0.408-0.295-0.603-0.461c-0.015-0.019-0.031-0.026-0.046-0.043l-0.665-0.664c-1.367-1.567-2.227-3.903-2.412-6.671C10.669,15.856,10.921,13.673,11.514,11.592')
        .attr({fill: "#000", stroke: "none"}).transform("t"+(width-30)+","+(height-30));  

    this.init = function(genomeId) {
        self.genomeId = genomeId;
        self.annotationsDisplayed.forEach(function(annotationDisplayed) {
            annotationDisplayed[1].remove();
        });
        self.annotationsDisplayed = [];
        for (var rank=0 ; rank < self.expressionLevels.length ; rank++)
            self.expressionLevels[rank].hide();
    }

    this.clear = function() {
        self.annotationsDisplayed.forEach(function(annotationDisplayed) {
            annotationDisplayed[1].remove();
        });
        self.annotationsDisplayed = [];
        for (var rank=0 ; rank < self.expressionLevels.length ; rank++)
            self.expressionLevels[rank].hide();
    }

    this.goUpDown = function(transY) {
       self.annotationsDisplayed.forEach(function(annotationDisplayed) {
            annotationDisplayed[1].transform("...t0,"+transY); 

            //we get the new y coord
            var transformedY = annotationDisplayed[1].attrs.y;

            annotationDisplayed[1].transform().forEach(function(t) {
               transformedY += t[2]; //we add all the translation values stored in the transform array linked to this object
            });

            if (transformedY < self.header) //important to not display the annotations mixed with the ticks
                annotationDisplayed[1].hide();
            else
                annotationDisplayed[1].show(); //don't forget to display back the annotations that have disappeared, if they're below the header space.  
        }); 
    };

    //this methods change the genomic range and check if some annotations displayed have to be removed or re-computed
    this.walkOnGenome = function(newGenomicCenter, annotation_id_to_blink) {
        for (var rank=0 ; rank < self.expressionLevels.length ; rank++)
            self.expressionLevels[rank].hide();
                
        if (annotation_id_to_blink)
            self.annotation_id_to_blink = annotation_id_to_blink;

        var genomicWidth = self.genomicRange[1] - self.genomicRange[0];
        self.genomicRange = [Math.round(newGenomicCenter-genomicWidth/2), Math.round(newGenomicCenter+genomicWidth/2)];

        if (self.genomicRange[0] <= 0) { //we cannot go below the chromosome start
            self.genomicRange[0] = 1;
            self.genomicRange[1] = width*self.genomicStep;
        }

        //we update the genomic sequence
        if (self.intervalBetweenTicks == 12) {
            self.displayGenomicSequence(self.current_genome.sequence.substring(self.genomicRange[0]-1, self.genomicRange[1]));
        }

        //we update the ticks
        self.ticks.forEach(function(tick) {
            tick.update();
        });

        //we reevaluate the annotationsDisplayed
        var current_ids =[], toBeRemoved = [];
        for (var i=0 ; i < self.annotationsDisplayed.length ; i++) {
            var annotation = self.annotationsDisplayed[i][0];
            if (annotation.genomicPositions[1] <= self.genomicRange[0] || annotation.genomicPositions[0] >= self.genomicRange[1])
                toBeRemoved.push(self.annotationsDisplayed[i]);
            else
              current_ids.push(annotation._id);
        }
        
        toBeRemoved.forEach(function(a) {
            a[1].remove();
            self.annotationsDisplayed.splice(self.annotationsDisplayed.indexOf(a),1);
        });
        
        //the remaining annotations are re-computed for their rect
        self.annotationsDisplayed.forEach(function(annotationDisplayed) {
            var startGenomicPos = annotationDisplayed[0].genomicPositions[0] >= self.genomicRange[0] ? annotationDisplayed[0].genomicPositions[0]: self.genomicRange[0]; 
            var length_sequence = annotationDisplayed[0].genomicPositions[1] - startGenomicPos+1,
                x = self.intervalBetweenTicks == 12 ? (startGenomicPos-self.genomicRange[0])*width/120 : (startGenomicPos-self.genomicRange[0])/self.genomicStep,
                y = annotationDisplayed[0].rank*25+self.header,
                w = self.intervalBetweenTicks == 12 ? length_sequence*width/120 : length_sequence/self.genomicStep,
                h = 20;

            if (w < 20) {
                annotationDisplayed[1].attr('path',"M"+x+","+y+"H"+(x+w)+"V"+(y+h)+"H"+(x)+"V"+y+"Z");   
            }
            else {
                if (annotationDisplayed[0].genomicStrand == '+') {
                    annotationDisplayed[1].attr('path', "M"+x+","+y+"H"+(x+w-10)+"L"+(x+w)+","+(y+h/2)+"L"+(x+w-10)+","+(y+h)+"H"+x+"V"+y+"Z");
                } else {
                    annotationDisplayed[1].attr('path', "M"+(x+10)+","+y+"H"+(x+w)+"V"+(y+h)+"H"+(x+10)+"L"+x+","+(y+h/2)+"L"+(x+10)+","+y+"Z"); 
                }
            }

            if (annotationDisplayed[0]._id == self.annotation_id_to_blink) {
                self.blinck(annotationDisplayed[0], annotationDisplayed[1]);
                self.annotation_id_to_blink = undefined;
            }
            //we get the new y coord
            var transformedY = annotationDisplayed[1].attrs.y;

            annotationDisplayed[1].transform().forEach(function(t) {
               transformedY += t[2]; //we add all the translation values stored in the transform array linked to this object
            });

            if (transformedY < self.header) //important to not display the annotations mixed with the ticks
                annotationDisplayed[1].hide();
            else
                annotationDisplayed[1].show();
        });

        //the browser asks the server for new annotations
        socket.send(JSON.stringify({
                'header':'genome browser dragged', 
                'genome_id': self.genomeId,
                'genomic_range': self.genomicRange,
                'current_ids': current_ids
            }));
    };

    this.getCurrentGenomicCenter = function() {
        var genomicWidth = self.genomicRange[1] - self.genomicRange[0];
        return self.genomicRange[0]+genomicWidth/2;
    };

    this.highlightGenomicSubSequence = function(start, end) {
        for (var i=1 ; i <= self.residues.length ; i++) {
            var residue = self.residues[i-1];
            if (i >= start && i <= end) {
                residue[1].attr({fill:'grey'});
            } else {
                if (residue[0].attr('text') == 'A')
                    residue[1].attr({fill:'#008000'});
                else if (residue[0].attr('text') == 'U' || residue[0].attr('text') == 'T')
                    residue[1].attr({fill:'#FFA500'});
                else if (residue[0].attr('text') == 'G')
                    residue[1].attr({fill:'#FF0000'});
                else if (residue[0].attr('text') == 'C')
                    residue[1].attr({fill:'#FF00FF'});
                else
                    residue[1].attr({fill:'black'});    
            }
        }
    };

    this.displayExpressionLevels = function (expressionLevels) {
        for (var rank=0 ; rank < self.expressionLevels.length ; rank++)
            if (!expressionLevels[rank]) //no expression at this rank
                self.expressionLevels[rank].update(0);
            else
                self.expressionLevels[rank].update(expressionLevels[rank]*100);
    };

    //display in the browser the annotation. This method creates the rectangle that will be used to display the genomic annotation.
    this.displayAnnotation = function(annotation) {
        var startGenomicPos = annotation.genomicPositions[0] >= self.genomicRange[0] ? annotation.genomicPositions[0]: self.genomicRange[0]; 
        var length_sequence = annotation.genomicPositions[1] - startGenomicPos+1;

        annotation.rank = 0

        for (var i=0 ; i < self.annotationsDisplayed.length ; i++) {
            var _annotation = self.annotationsDisplayed[i][0];
            if (_annotation._id == annotation._id) //just to be sure that it is not already displayed
                return;
            if (annotation.genomicPositions[0] >=  _annotation.genomicPositions[0] && annotation.genomicPositions[0] <=  _annotation.genomicPositions[1] ||
                annotation.genomicPositions[1] >=  _annotation.genomicPositions[0] && annotation.genomicPositions[1] <=  _annotation.genomicPositions[1] ||
                _annotation.genomicPositions[0] >=  annotation.genomicPositions[0] && _annotation.genomicPositions[0] <=  annotation.genomicPositions[1] ||
                _annotation.genomicPositions[1] >=  annotation.genomicPositions[0] && _annotation.genomicPositions[1] <=  annotation.genomicPositions[1] ) {
                if (annotation.rank <= _annotation.rank)
                    annotation.rank = _annotation.rank+1;
            }    
        }

        var rect = undefined,
            x = self.intervalBetweenTicks == 12 ? (startGenomicPos-self.genomicRange[0])*width/120 : (startGenomicPos-self.genomicRange[0])/self.genomicStep,
            y = annotation.rank*25+self.header,
            w = self.intervalBetweenTicks == 12 ? length_sequence*width/120 : length_sequence/self.genomicStep,
            h = 20;

        if (w < 20)
            rect = self.canvas.path("M"+x+","+y+"H"+(x+w)+"V"+(y+h)+"H"+(x)+"V"+y+"Z");   
        else {
            if (annotation.genomicStrand == '+') {
                rect = self.canvas.path("M"+x+","+y+"H"+(x+w-10)+"L"+(x+w)+","+(y+h/2)+"L"+(x+w-10)+","+(y+h)+"H"+x+"V"+y+"Z");
            } else {
                rect = self.canvas.path("M"+(x+10)+","+y+"H"+(x+w)+"V"+(y+h)+"H"+(x+10)+"L"+x+","+(y+h/2)+"L"+(x+10)+","+y+"Z");
            }
        }
        
        self.annotationsDisplayed.push([annotation, rect]);               

        rect.click(function() {
            var startGenomicPos = annotation.genomicPositions[0] >= self.genomicRange[0] ? annotation.genomicPositions[0]: self.genomicRange[0];
            var length_sequence = annotation.genomicPositions[1] - startGenomicPos+1; 
            self.highlightGenomicSubSequence(startGenomicPos-self.genomicRange[0]+1, startGenomicPos-self.genomicRange[0]+length_sequence);
            socket.send(JSON.stringify({
                'header':'get annotation to highlight', 
                'annotation': annotation
            }));
        });


        rect.attr({
            title : "Product: "+annotation.product +"\nClass: "+ annotation.class  +"\nLength: "+ (annotation.genomicPositions[1]-annotation.genomicPositions[0]+1) +"nts\nStrand: " +   
            annotation.genomicStrand + "\nOrganism: "+ annotation.organism       
        }).toBack(); //toBack() is to have the genomic annotations displayed below some widgets like the dock with the icons
        
        var color = self.getAnnotationColor(annotation);
        
        if (annotation._id == self.annotation_id_to_blink) {
            rect.attr({
                fill: color,                
                "fill-opacity": .4,
                stroke:color
            });

            self.blinck(annotation, rect);
            self.annotation_id_to_blink = undefined;

        } else {
            rect.attr({
                fill:'white',
                "fill-opacity": .4,
                stroke:'white'
            });    
            rect.animate({fill: color , stroke: color}, 1000);
        }

    };

    this.scale = function(ratio) {
        if (ratio <1 && self.intervalBetweenTicks <= 12) // no minus zoom, otherwise we will have less than 1 nt between two consecutive ticks
            return;

        for (var rank=0 ; rank < self.expressionLevels.length ; rank++)
            self.expressionLevels[rank].hide();

        var scale = (self.genomicRange[1] - self.genomicRange[0])*ratio;

        var genomicCenter = self.getCurrentGenomicCenter();

        self.genomicRange[0] = genomicCenter-scale/2;
        self.genomicRange[1] = genomicCenter+scale/2;

        if (self.genomicRange[0] <= 0) { //we cannot go below the chromosome start
            self.genomicRange[0] = 1;
            self.genomicRange[1] = scale;
        }

        self.intervalBetweenTicks = Math.round((self.genomicRange[1]-self.genomicRange[0])/10);

        if (self.intervalBetweenTicks <= 12) { //we cannot magnify plus below 12nts between each ticks
            self.intervalBetweenTicks = 12;            
            self.genomicRange[1] = self.genomicRange[0]+120; //we want 120 nts between ends (12nts*10 ticks)
            self.residues.forEach(function(residue) {
                residue[0].attr({opacity:1});
                residue[1].attr({opacity:1});
            });        
        } else {
            self.residues.forEach(function(residue) {
                residue[0].attr({opacity:0});
                residue[1].attr({opacity:0});
            });     
        }
        
        self.genomicStep = (self.genomicRange[1] - self.genomicRange[0]+1)/width;

        self.walkOnGenome(self.getCurrentGenomicCenter());
       
    };

    this.displayGenomicSequence = function(genomicSequence) {
        if (self.intervalBetweenTicks == 12) {
            //console.log(genomicSequence);
            var chars = genomicSequence.split('');
            for (var i = 0 ; i < chars.length-1 ; i++) {
                var residue_char = chars[i];
                self.residues[i][0].attr({'text':residue_char});
                if (residue_char == 'A')
                    self.residues[i][1].attr({fill:'#008000'});
                else if (residue_char == 'U' || residue_char == 'T')
                    self.residues[i][1].attr({fill:'#FFA500'});
                else if (residue_char == 'G')
                    self.residues[i][1].attr({fill:'#FF0000'});
                else if (residue_char == 'C')
                    self.residues[i][1].attr({fill:'#FF00FF'});
                else
                    self.residues[i][1].attr({fill:'black'});    
            }
        }
    };

    this.annotationValidated = function(annotation) {
        for (var i=0 ; i < self.annotationsDisplayed.length ; i++) {
            var _annotation = self.annotationsDisplayed[i][0];
            if (_annotation._id == annotation._id) {
                self.annotationsDisplayed[i][1].attr({
                    fill: 'orange',                
                    stroke:'orange'
                });
                break;
            }  
        }    
    };

    this.annotationInValidated = function(annotation) {
        for (var i=0 ; i < self.annotationsDisplayed.length ; i++) {
            var _annotation = self.annotationsDisplayed[i][0];
            if (_annotation._id == annotation._id) {
                self.annotationsDisplayed[i][1].attr({
                    fill: self.getAnnotationColor(annotation),                
                    stroke: self.getAnnotationColor(annotation)
                });
                break;
            }           
        }    
    };

    this.blinck = function(annotation, rect) {
        rect.animate({fill:'white', stroke:'white'}, 500, 
            function() {
                rect.animate({fill:self.getAnnotationColor(annotation), stroke:self.getAnnotationColor(annotation)}, 500, function() {
                    rect.animate({fill:'white', stroke:'white'}, 500, function() {
                        rect.animate({fill:self.getAnnotationColor(annotation), stroke:self.getAnnotationColor(annotation)}, 500);
                    });
                });
            }
        );
    };

    this.getAnnotationColor = function(annotation) {
        if (annotation.class.toUpperCase() == "TRNA")
            return 'purple';
        else if (annotation.class.toUpperCase() == "RRNA")
            return 'red';
        else if (annotation.class.toUpperCase() == "INTRON")
            return 'orange';
        else if (annotation.class.toUpperCase() == "GENE")
            return color = '#0B6121';
        else if (annotation.class.toUpperCase() == "MRNA")
            return color = '#088A08';
        else if (annotation.class.toUpperCase() == "CDS")
            return color = '#A5DF00';
        else if (annotation.class.toUpperCase() == "NCRNA")
            return color = '#0040FF'; 
        else 
            return color = 'grey';     
    };

};