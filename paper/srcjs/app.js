import CanvasImage from "./canvasImage";
import React, { useState, useEffect } from 'react';
import "./app.css";

const App = (props) => {

  console.log(props)

  return (
    <>
      <CanvasImage
        data={props.data}
        width={props.width}
        height={props.height}
      />
    </>
  )
};

export default App;
