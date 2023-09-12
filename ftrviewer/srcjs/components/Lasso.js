import React, { useRef, useEffect } from 'react'
import Paper from 'paper'
import Draw from './Draw'

const Lasso = (props) => {
  const canvasRef = useRef(null)

  useEffect(() => {
    const canvas = canvasRef.current
    Paper.setup(canvas)
    Draw(props)
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  return (
    <canvas
      ref={canvasRef}
      id="canvas"
      className="overlay"
      resize="true"
      style={{
        height: props.height,
        width: props.width,
        position: 'absolute',
        top: 0,
        left: "15px",
      }}
    />
  )
}

export default Lasso
