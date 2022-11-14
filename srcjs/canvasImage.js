import React, { useRef, useEffect, useState } from 'react';
import Paper from 'paper';
import { Key, Point } from 'paper/dist/paper-core';
import "./canvasImage.css";

window.MyLib = {};

const CanvasImage = (props) => {

  var rawdata = [];
  for (const key of Object.keys(props.data)) {
    rawdata.push(props.data[key]);
  }

  // fetch images
  var images = rawdata.map((im, i) => {
    var arr = new Uint8ClampedArray(im.data);
    let imageData = new ImageData(arr, im.dimx, im.dimy);
    return imageData;
  });

  // Create state for active images
  const [active, setActive] = useState(images.map(() => true));

  // items
  var items = images.map((im, i) => {
    return "heimage" + i;
  });

  // Create a reference to the paper js canvas element
  const canvasRefImage = useRef(null);
  // Setup paper js canvas
  useEffect(() => {
    const canvas = canvasRefImage.current;
    Paper.setup(canvas);
  }, []);

  // Create a reference to the raster image elements
  const canvasRefs = useRef([]);
  canvasRefs.current = [];

  // Put image in context
  const draw = (ctx, imageData, i) => {
    //var arr = new Uint8ClampedArray(props.data.data);
    //let imageData = new ImageData(arr, props.data.dimx, props.data.dimy);
    //var imageData = images[i];
    ctx.putImageData(imageData, 0, 0);
  };

  // Add ellement to the canvasRefs array
  const addToRefs = (el) => {
    if (el && !canvasRefs.current.includes(el)) {
      canvasRefs.current.push(el);
    }
  };

  // draw images in context
  useEffect(() => {
    images.map((im, i) => {
        const canvas = canvasRefs.current[i];
        const ctx = canvas.getContext('2d');
        canvas.width = im.width;
        canvas.height = im.height;
        draw(ctx, im, i);
    });
  }, [draw]);


  // Place images in paper js canvas
  useEffect(() => {

    // Read rasters from DOM img elements
    var rasters = items.map((item, index) => {
        MyLib[index] = {id: item, index: index, angle: 0, scalefactor: 1, shift_x: 0, shift_y: 0, selected: false, visible: true, flip: false};
        var domimg = document.getElementById(item);
        let raster = new Paper.Raster(domimg);
        raster.onMouseEnter = () => {
            Paper.view.element.style.setProperty("cursor", "move");
        };
        raster.onMouseLeave = () => {
            Paper.view.element.style.setProperty("cursor", null);
        };
        raster.position = Paper.view.center;
        raster.scalefactor = 1;
        return raster;
    });

    // Create rectangle for HE image borders
    var rectBorder = new Paper.Path.Rectangle({
        point: rasters[0].bounds.topLeft,
        size: rasters[0].bounds.size,
        strokeColor: 'black',
        fillColor: 'white'
    });
    rectBorder.dashArray = [10, 4];
    rectBorder.sendToBack();

    // Create rectangle for HE image
    const paddPoint = new Paper.Point(20, 20);
    const newTopLeftCorner = rasters[0].bounds.topLeft.subtract(paddPoint);
    const newSize = rasters[0].bounds.size.add(paddPoint.multiply(2));
    var rect = new Paper.Path.Rectangle({
        point: newTopLeftCorner,
        size: newSize,
        strokeColor: 'white',
        fillColor: "white"
    });
    rect.sendToBack();

    // Create background
    var background = new Paper.Path.Rectangle({
        point: [0, 0],
        size: Paper.view.size,
        strokeColor: "#E6E6E6",
        fillColor: "#F0F0F0"
    });
    background.sendToBack();

    // Add event handlers to rasters
    window.globals = [];
    rasters.map((raster, index) => {
        if (raster.loaded) {

            var corners = {
                topRight: raster.bounds.topRight.subtract(raster.position),
            };

            var diffAngle = corners.topRight.getDirectedAngle(raster.bounds.rightCenter.subtract(raster.position));
            var isFlipped = false;

            var topRight = new Paper.Path.Circle({
                center: raster.position.add(corners.topRight),
                radius: 5,
                fillColor: "#4682B4",
                strokeColor: "black",
                onMouseEnter: function () {
                  Paper.view.element.style.setProperty("cursor", "grab");
                },
                onMouseLeave: function () {
                  Paper.view.element.style.setProperty("cursor", null);
                },
              });

            raster.onMouseDown = (e) => {
                raster.bringToFront();
                topRight.bringToFront();
                if (Key.isDown("q")) {
                    raster.opacity = Math.max(raster.opacity - 0.15, 0.1);
                }
                if (Key.isDown("w")) {
                    raster.opacity = Math.min(raster.opacity + 0.15, 1);
                }
                if (Key.isDown("r")) {
                  raster.position = Paper.view.center;
                  raster.rotation = 0;
                  raster.scale(1/raster.scalefactor);
                  raster.scalefactor = 1;
                  raster.opacity = 1;
                  topRight.position = raster.bounds.topRight;
                }
                if (Key.isDown("f")) {
                  raster.scale(-1, 1);
                  isFlipped = !isFlipped;
                  var diffX = (topRight.position.x - raster.position.x);
                  topRight.position.x = raster.position.x - diffX;
                  window.MyLib[index].flip = !window.MyLib[index].flip;
                }
                raster.diff = e.point.subtract(raster.position);
            };

            raster.onMouseDrag = (e) => {
                raster.position = e.point.subtract(raster.diff);
                var curAngle = raster.rotation - (isFlipped ? -diffAngle : diffAngle);
                var newTopRight = new Paper.Point({
                    angle: curAngle,
                    length: corners.topRight.length*raster.scalefactor
                })
                topRight.position = raster.position.add(newTopRight);
            };

            raster.onMouseUp = (e) => {
              window.MyLib[index].angle = raster.rotation;
              window.MyLib[index].scalefactor = raster.scalefactor;
              var shift_along_xy = Paper.view.center.subtract(raster.position);
              window.MyLib[index].shift_x = shift_along_xy.x;
              window.MyLib[index].shift_y = shift_along_xy.y;
              Object.keys(window.MyLib).forEach(key => {
                if (key != index) {
                  window.MyLib[key].selected = false;
                }
              });
              window.MyLib[index].selected = true;
            };

            topRight.onMouseDrag = (e) => {
                if (Key.isDown("shift")) {
                    var curPoint = e.point.subtract(raster.position);
                    var prevPoint = topRight.position.subtract(raster.position);
                    var angle = -curPoint.getDirectedAngle(prevPoint);
                    raster.rotate(angle);
                }
                else {
                    var curDistance = raster.position.getDistance(e.point);
                    var prevDistance = corners.topRight.rotate(raster.rotation).length;
                    var scalefactor = curDistance/prevDistance;
                    raster.scale(scalefactor/raster.scalefactor);
                    raster.scalefactor = scalefactor;
                }

                var curAngle = raster.rotation - (isFlipped ? -diffAngle : diffAngle);
                var newTopRight = new Paper.Point({
                    angle: curAngle,
                    length: corners.topRight.length*raster.scalefactor
                })
                topRight.position = raster.position.add(newTopRight);

            };

            topRight.onMouseUp = (e) => {
              window.MyLib[index].angle = raster.rotation;
              window.MyLib[index].scalefactor = raster.scalefactor;
              var shift_along_xy = Paper.view.center.subtract(raster.position);
              window.MyLib[index].shift_x = shift_along_xy.x;
              window.MyLib[index].shift_y = shift_along_xy.y;
              Object.keys(window.MyLib).forEach(key => {
                if (key != index) {
                  window.MyLib[key].selected = false;
                }
              });
              window.MyLib[index].selected = true;
            };

            window.globals[index] = {
              visible: true,
              hide: () => {
                raster.visible = !raster.visible;
                topRight.visible = !topRight.visible;
              }
            }

        }
    });

  }, []);

  return (
    <>
      <div className="image-buttons">
        {images.map((image, i) => {
          return (
            <button
              className="button-4"
              key={i}
              style={{
                backgroundColor: active[i] ? "#dbf1f4" : "#E6E6E6",
                color: active[i] ? "black" : "#C8C8C8",
              }}
              onClick={() => {
                window.globals[i].hide();
                const copy = [...active];
                copy[i] = !copy[i];
                setActive(copy);
              }}
            >
              Image {i}
            </button>
          );
        })}
      </div>
      <div className="canvas-container">
        <canvas
            ref={canvasRefImage}
            id="canvas-annotate"
            style={{
                height: "600px",
                width: "100vh",

            }}
        />
        {
          items.map((item, index) => {
              return (
                  <canvas
                      ref={addToRefs}
                      key={index}
                      style={{ display: "none" }}
                      id={"heimage" + index}
                  ></canvas>
              );
          })
        }
      </div>
    </>
  );
}

export default CanvasImage;
