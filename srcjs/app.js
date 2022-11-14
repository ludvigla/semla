import CanvasImage from "./canvasImage";
import React, { useState, useEffect } from 'react';
import "./app.css";

const App = (props) => {

  return (
    <>
      <CanvasImage
        data={props.data}
      />
    </>
  )
};

export default App;
