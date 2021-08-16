class DirBase {
  constructor(x, y) {
    this.x = x;
    this.y = y;
  }

  x = undefined;
  y = undefined;
}

class Direction extends DirBase {
  get() { return Object.assign(Object.create(Object.getPrototypeOf(this)), this); }
  getInv() {
    var temp = this.get();
    
    temp.x *= -1;
    temp.y *= -1;

    return temp;
  }
}

var diffDirection = (dx, dy) => {
  var _dx = parseFloat(dx), _dy = parseFloat(dy);
  if ([_dx, _dy].every(v => v == 0) || [_dx, _dy].some(n => Number.isNaN(n))) {
    throw new Error(`bad input; dx: \`${dx}\`, dy: \`${dy}\``);
  }
  
  var d = Math.sqrt(_dx ** 2 + _dy ** 2);
  return new Direction(dx / d, dy / d);
}

var polarDirection = theta => new Direction(Math.cos(theta), Math.sin(theta));

class Member {
  constructor(src, dst, direction) {
    this.src = src;
    this.dst = dst;
    this.direction = direction;
  }

  src = undefined;
  dst = undefined;
  direction = undefined;
}

class Force extends DirBase {}

var polarForce = (r, theta) => new Force(r * Math.cos(theta), r * Math.sin(theta));
var xForce = f => new Force(f, 0), yForce = f => new Force(0, f);

class Reaction {
  constructor(name, jointName, direction) {
    this.name = name;
    this.jointName = jointName;
    this.direction = direction;
  }

  name = undefined;
  jointName = undefined;
  direction = undefined;
}


var gauss = A => {
  var n = A.length;

  for (var i = 0; i < n; ++i) {
    // Search for maximum in this column
    var maxEl = Math.abs(A[i][i]);
    var maxRow = i;
    for (var k = i + 1; k < n; ++k) {
      if (Math.abs(A[k][i]) > maxEl) {
        maxEl = Math.abs(A[k][i]);
        maxRow = k;
      }
    }

    // Swap maximum row with current row (column by column)
    for (var k = i; k < n+ 1; ++k) {
        var tmp = A[maxRow][k];
        A[maxRow][k] = A[i][k];
        A[i][k] = tmp;
    }

    // Make all rows below this one 0 in current column
    for (k = i + 1; k < n; ++k) {
      var c = -A[k][i] / A[i][i];
      for (var j = i; j < n + 1; ++j) {
        if (i == j) {
          A[k][j] = 0;
        } else {
          A[k][j] += c * A[i][j];
        }
      }
    }
  }

  // Solve equation Ax=b for an upper triangular matrix A
  var x = new Array(n);
  for (var i = n - 1; i > -1; --i) {
    x[i] = A[i][n] / A[i][i];

    for (var k = i - 1; k > -1; --k) {
      A[k][n] -= A[k][i] * x[i];
    }
  }

  return x;
}

var solve = (memberList, forceList, reaction1, reaction2, reaction3) => {
	var jointRelation = {}, nMembers = memberList.length;

  // load in member list
  var memberNames = [];
	for (var member of memberList) {
    member.src = member.src.replace(/-/g, "_");
    member.dst = member.dst.replace(/-/g, "_");

		if (jointRelation[member.src] === undefined) jointRelation[member.src] = {};
    if (jointRelation[member.dst] === undefined) jointRelation[member.dst] = {};

    if (member.src > member.dst) {
      member.direction = member.direction.getInv();
      [member.src, member.dst] = [member.dst, member.src];
    }

    jointRelation[member.src][member.dst] = member.direction.get();
    jointRelation[member.dst][member.src] = member.direction.getInv();

    memberNames.push(`${member.src}-${member.dst}`);
	}
  memberNames = memberNames.concat([reaction1, reaction2, reaction3].map(v => v.name)).sort();

  if (memberNames.length !== Object.keys(jointRelation).length * 2) {
    throw new Error("wrong number of joints/members/reactions");
  }

  // construct simultaneous equations
  var equations = {};
  for (var jointName in jointRelation) {
    equations[jointName] = {
      x: memberNames.reduce((r, v) => { r.entries[v] = 0; return r; }, { entries: {}, sum: 0 }),
      y: memberNames.reduce((r, v) => { r.entries[v] = 0; return r; }, { entries: {}, sum: 0 }),
    };

    // add the members to the equation(s)
    for (var k in jointRelation[jointName]) {
      var mod = [jointName, k].sort().join("-");
      equations[jointName].x.entries[mod] = jointRelation[jointName][k].x;
      equations[jointName].y.entries[mod] = jointRelation[jointName][k].y;
    }
  }

  // add forces to simultaneous equations
  for (var [jointName, force] of forceList) {
    // `sum` is in the opposite direction to the force => `-=`
    equations[jointName].x.sum -= force.x;
    equations[jointName].y.sum -= force.y;
  }

  // add reaction forces
  for (var reaction of [reaction1, reaction2, reaction3]) {
    equations[reaction.jointName].x.entries[reaction.name] = reaction.direction.x;
    equations[reaction.jointName].y.entries[reaction.name] = reaction.direction.y;
  }

  // solve equations
  var entrySortFunc = e => Object.entries(e).sort(
    ([k1], [k2]) => memberNames.indexOf(k1) - memberNames.indexOf(k2)
  ).map(v => v[1]);
  var matrix = Object.values(equations).reduce((r, v) => ([
    ...r,
    [...entrySortFunc(v.x.entries), v.x.sum],
    [...entrySortFunc(v.y.entries), v.y.sum],
  ]), []);
  var res = gauss(matrix);

	return Object.fromEntries(res.map((_, i) => [memberNames[i], res[i]]));
};
