const express = require("express");
const bodyParser = require("body-parser");
const cors = require("cors");

const calculator = require("./build/Release/my_addon");

const app = express();
app.use(cors());
app.use(bodyParser.urlencoded({ extended: false }));
app.use(bodyParser.json());

const port = 3000;

app.post("/calculate", (request, response) => {
  const {
    points,
    A = 40,
    B = 45,
    C = 10,
    DD = 45,
    omega = 55,
    k = 50,
    V = 40,
    S = 1300,
    L = 1.33,
    c = 0.28,
    N = 21000,
    h = 5,
    sigma = -0.5,
    a = 0.1,
    b = 0.1,
    p = 0.1,
    alpha = 6.611,
    betta = 3.636,
    x = 0,
    y = 1,
    koeffXMin = 50,
    koeffXMax = 50,
    koeffYMin = 50,
    koeffYMax = 50
  } = request.body;

  try {
    const values = calculator.calculate(
      points,
      A,
      B,
      C,
      DD,
      omega,
      k,
      V,
      S,
      L,
      c,
      N,
      h,
      sigma,
      a,
      b,
      p,
      alpha,
      betta,
      x,
      y,
      koeffXMin,
      koeffXMax,
      koeffYMin,
      koeffYMax
    );

    return response
      .json({ values })
      .status(200)
      .end();
  } catch (error) {
    console.log(error);
    return response
      .json(error)
      .status(500)
      .end();
  }
});

app.listen(port, () => {
  console.log(`server started on http://localhost:${port}`);
});