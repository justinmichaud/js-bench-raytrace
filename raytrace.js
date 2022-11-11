let rnd = null
let useDelete = false

class Ray {
    constructor(pos, dir) {
        this.pos = pos
        this.dir = dir
    }

    transform(t, s) {
        return new Ray(add(scale(this.pos, s), t), scale(this.dir, s))
    }
}

class Point {
    constructor(x, y, z) {
        this.x = x
        this.y = y
        this.z = z
        this.a = 0
        if (useDelete)
            delete this.a
    }

    get r() { return this.x }
    get g() { return this.y }
    get b() { return this.z }
}

class Intersection {
    constructor(pos, normal) {
        this.pos = pos
        this.normal = normal
    }

    transform(t, s) {
        let pos = add(scale(this.pos, s), t)
        // This must be a uniform scale, since the normal is a covector
        let normal = scale(this.normal, s)
        return new Intersection(pos, normal)
    }
}

class ConstantMaterial {
    constructor(reflectance, emittance) {
        this.reflectance = reflectance
        this.emittance = emittance
    }
}

class Sphere {
    constructor(r, pos) {
        this.r = r
        this.pos = pos
    }

    intersect(ray) {
        ray = ray.transform(scale(this.pos, -1.0/this.r), 1.0/this.r)
        let {pos, dir} = ray

        let A = dot(dir,dir)
        let B = 2*dot(dir,pos)
        let C = dot(pos,pos)-1

        let ts = quadraticRoots(A,B,C)
        if (ts.length == 0)
            return null
        
        let t = ts[0]
        if (t < 0 || (ts[1] >= 0 && ts[1] < t)) t = ts[1]
        if (t < 0)
            return null
        let modelSpace = new Intersection(add(pos, scale(dir,t)), add(pos, scale(dir,t)))
        return modelSpace.transform(this.pos, this.r)
    }
}

class Cube {
    constructor(r, pos, invertNormals = false) {
        this.r = r
        this.pos = pos
        this.invertNormals = invertNormals
    }

    intersect(ray) {
        ray = ray.transform(scale(this.pos, -1.0/this.r), 1.0/this.r)
        let {pos, dir} = ray

        const faces = [
            [new Point(1,0,0), new Point(0,1,0), new Point(1,0,0)],
            [new Point(0,0,0), new Point(0,1,0), new Point(-1,0,0)],
            [new Point(0,1,0), new Point(1,0,0), new Point(0,1,0)],
            [new Point(0,0,0), new Point(1,0,0), new Point(0,-1,0)],
            [new Point(0,0,1), new Point(1,0,0), new Point(0,0,1)],
            [new Point(0,0,0), new Point(1,0,0), new Point(0,0,-1)],
        ]

        let closest = null

        for (let [point, px, norm] of faces) {
            if (dot(norm, dir) == 0)
                continue
            let t = dot(add(point, scale(pos, -1)),norm)/dot(norm,dir)
            if (t < 0)
                continue
            let p = add(pos, scale(dir,t))
            let py = cross(px, norm)
            py.x = Math.abs(py.x)
            py.y = Math.abs(py.y)
            py.z = Math.abs(py.z)
            let planeCoords = [dot(p,px), dot(p,py)]
            if (planeCoords[0] < 0 || planeCoords[0] > 1 || planeCoords[1] < 0 || planeCoords[1] > 1)
                continue
            if (!closest || t < closest.t) {
                closest = {p, norm, t}
            }
        }

        if (!closest)
            return null
        let {p, norm} = closest
        if (this.invertNormals)
            norm = scale(norm, -1)
        let modelSpace = new Intersection(p, norm)
        return modelSpace.transform(this.pos, this.r)
    }
}

class WorldObject {
    constructor(shape, material) {
        this.shape = shape
        this.material = material
    }
}

function quadraticRoots(A, B, C) {
    // Computes the Citardauq Formula
    let D = 0.0
    let q = 0.0
    let roots = []

	if( A == 0 ) {
		if( B == 0 )
			return []
		else
			return [-C/B]
	} else {
		// Compute the discriminant D=b^2 - 4ac
        D = B*B - 4*A*C
		if( D < 0 )
			return []
		else {
            let signB = Math.sign(B)
            if (B == 0)
                signB = 1.0
			q = -(B + signB*Math.sqrt(D)) / 2.0
            roots[0] = q / A
			if( q != 0 ) {
				roots[1] = C / q
			}
			return roots
		}
	}
}

function dot(a, b) {
    return a.x*b.x+a.y*b.y+a.z*b.z
}

function cross(a, b) {
    let r = new Point(a.x,a.y,a.z)
    r.x = a.y*b.z-a.z*b.y
    r.y = a.z*b.x-a.x*b.z
    r.z = a.x*b.y-a.y*b.x
    return r
}

function scale(a, s) {
    a = new Point(a.x,a.y,a.z)
    a.x *= s
    a.y *= s
    a.z *= s
    return a
}

function magnitude(r) {
    return Math.sqrt(dot(r,r))
}

function distance(a, b) {
    return magnitude(add(a, scale(b, -1.0)))
}

function normalize(r) {
    return scale(r, 1.0 / magnitude(r))
}

function add(a, b) {
    a = new Point(a.x,a.y,a.z)
    a.x += b.x
    a.y += b.y
    a.z += b.z
    return a
}

function sub(a, b) {
    a = new Point(a.x,a.y,a.z)
    a.x -= b.x
    a.y -= b.y
    a.z -= b.z
    return a
}

// This produces a random sample from a cosine-weighted hemisphere with the given normal.
// https://www.scratchapixel.com/lessons/3d-basic-rendering/global-illumination-path-tracing/global-illumination-path-tracing-practical-implementation
function sampleHemisphere(normal) {
    let r1 = rnd()
    let r2 = rnd()

    let theta = Math.acos(Math.sqrt(r1))
    let sinTheta = Math.sin(theta)
    let phi = 2 * Math.PI * r2
    let local = new Point(sinTheta*Math.cos(phi), Math.cos(theta), sinTheta*Math.sin(phi))
                             
    let hemiSpaceX = null
    if (Math.abs(normal.x) > Math.abs(normal.y))
        hemiSpaceX = normalize(new Point(normal.z, 0, -normal.x))
    else
        hemiSpaceX = normalize(new Point(0, -normal.z, normal.y))
    let hemiSpaceY = normal
    let hemiSpaceZ = cross(hemiSpaceX, hemiSpaceY)

    return normalize(add(add(scale(hemiSpaceX, local.x), scale(hemiSpaceY, local.y)), scale(hemiSpaceZ, local.z)))
}

const BRIGHT = 100.0
const GREEN = new ConstantMaterial(new Point(0.13, 0.54, 0.13), new Point(0,0,0))
const WHITE = new ConstantMaterial(new Point(1,1,1), new Point(0,0,0))
const LIGHT = new ConstantMaterial(new Point(0, 0, 0), new Point(0.5,0.5,0.5))

const WORLD = [
    new WorldObject(new Sphere(2.5, new Point(-4,-7.5,-13)), GREEN),
    new WorldObject(new Cube(4, new Point(2,-8,-13)), GREEN),
    new WorldObject(new Cube(20, new Point(-10,-10,-19), true), WHITE),
    new WorldObject(new Cube(10, new Point(-5,8,-12)), LIGHT),
]

function intersectWorld(ray) {
    let closestIntersection = null

    for (let o of WORLD) {
        let i = o.shape.intersect(ray)
        if (!i)
            continue
        let { pos, normal } = i
        normal = normalize(normal)

        if (!closestIntersection || distance(closestIntersection.pos, ray.pos) > distance(pos, ray.pos))
            closestIntersection = {pos, normal, o}
    }

    return closestIntersection
}

function colorForRay(ray, depth = 0) {
    let i = intersectWorld(ray)
    if (!i || depth >= 4)
        return new Point(0,0,0)
    ++depth
    let { pos, normal, o } = i

    const INNER = 1
    let result = new Point(0,0,0)

    for (let c=0; c<INNER; ++c) {
        let newDir = sampleHemisphere(normal)
        let newRay = new Ray(add(pos, scale(newDir, 0.1)), newDir)
        let cos_theta = dot(newRay.dir, normal)
        let incoming = colorForRay(newRay, depth)
        incoming = scale(incoming, cos_theta)
        result = add(result, incoming)
    }

    let BRDF = scale(o.material.reflectance, 1/INNER)
    result.x *= BRDF.x
    result.y *= BRDF.y
    result.z *= BRDF.z

    return add(scale(o.material.emittance, BRIGHT), result)
}

function makeRay(x, y, ar) {
    const fovy = 90
    let dir = new Point(x*ar+(2*rnd()-1)*0.001, -y+(2*rnd()-1)*0.001, -1.0/Math.tan(fovy*3.14159/180.0/2.0))
    return new Ray(new Point(0,0,0), normalize(dir))
}

export function setRnd(_rnd) {
    rnd = _rnd
}

export function setDelete(_useDelete) {
    useDelete = _useDelete
}

function trace(x, y) {
    let res = new Point(0,0,0)
    let INNER = 1
    for (let i=0; i<INNER; ++i)
        res = add(res, colorForRay(makeRay(x,y,1.0)))
    const gamma = 1.0 / 2.2
    res.x = Math.pow(res.x / INNER, gamma)
    res.y = Math.pow(res.y / INNER, gamma)
    res.z = Math.pow(res.z / INNER, gamma)
    return res
}

export function* generateImage(canvas) {
    let ctx = canvas.getContext('2d')
    let width = canvas.width
    let height = canvas.height
    let imageData = ctx.getImageData(0, 0, width, height)

    // https://hacks.mozilla.org/2011/12/faster-canvas-pixel-manipulation-with-typed-arrays/
    let raw = new Uint32Array(width*height*3)
    let output = new Uint8Array(imageData.data.buffer)
    let samples = 1

    while (true) {
        for (let y = 0; y < height; ++y) {
            for (let x = 0; x < width; ++x) {
                let i = (y * width + x)*3
                let i2 = (y * width + x)*4
                let {r, g, b} = trace(x/width*2.0-1.0, y/height*2.0-1.0)

                output[i2] = Math.min(255*((raw[i] = raw[i] + r) / samples), 255)
                output[i2+1] = Math.min(255*((raw[i+1] = raw[i+1] + g) / samples), 255)
                output[i2+2] = Math.min(255*((raw[i+2] = raw[i+2] + b) / samples), 255)
                output[i2+3] = 255
            }
        }

        samples++
        ctx.putImageData(imageData, 0, 0)
        yield true
    }
}