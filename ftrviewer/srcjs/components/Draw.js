import Paper from "paper";

const Draw = (props) => {

console.log(props)

  var myPath;
  var myPathErase;

  Paper.view.onMouseDown = (event) => {
    if(Paper.Key.isDown('shift')) {
      myPathErase = new Paper.Path();
      myPathErase.strokeColor = "#771122";
      myPathErase.fillColor = "#771122";
      myPathErase.strokeWidth = 5;
      myPathErase.opacity = 0.5;
    } else {
      myPath = new Paper.Path();
      myPath.strokeColor = "#117777";
      myPath.fillColor = "#117777";
      myPath.strokewidth = 5;
      myPath.opacity = 0.5;
    }
  };

  Paper.view.onMouseDrag = (event) => {
    if(Paper.Key.isDown('shift')) { 
      myPathErase.add(event.point)

    } else {
      myPath.add(event.point);
    }
  };

  Paper.view.onMouseUp = (event) => {

    if(Paper.Key.isDown('shift')) { 
      myPathErase.closed = true;
      props.passEraseData(myPathErase);
      myPathErase.remove();
    } else {
      myPath.closed = true;
      props.passData(myPath);
      myPath.remove()
    }
  };

  Paper.view.draw();
};

export default Draw;