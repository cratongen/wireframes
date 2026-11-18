/// Closed axis-aligned box. The empty box has `min > max`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Aabb {
    pub min: [f64; 3],
    pub max: [f64; 3],
}

impl Aabb {
    pub fn empty() -> Aabb {
        Aabb {
            min: [f64::INFINITY; 3],
            max: [f64::NEG_INFINITY; 3],
        }
    }

    pub fn from_points(points: &[[f64; 3]]) -> Aabb {
        points
            .iter()
            .fold(Aabb::empty(), |b, &p| b.union(&Aabb { min: p, max: p }))
    }

    pub fn union(&self, o: &Aabb) -> Aabb {
        Aabb {
            min: [0, 1, 2].map(|k| self.min[k].min(o.min[k])),
            max: [0, 1, 2].map(|k| self.max[k].max(o.max[k])),
        }
    }

    pub fn intersects(&self, o: &Aabb) -> bool {
        (0..3).all(|k| self.min[k] <= o.max[k] && o.min[k] <= self.max[k])
    }

    pub fn contains(&self, p: [f64; 3]) -> bool {
        (0..3).all(|k| self.min[k] <= p[k] && p[k] <= self.max[k])
    }

    pub fn center(&self) -> [f64; 3] {
        [0, 1, 2].map(|k| 0.5 * (self.min[k] + self.max[k]))
    }

    /// Squared distance from `p` to the box, 0 inside.
    pub fn distance2(&self, p: [f64; 3]) -> f64 {
        (0..3)
            .map(|k| {
                let d = (self.min[k] - p[k]).max(p[k] - self.max[k]).max(0.0);
                d * d
            })
            .sum()
    }

    pub fn expand(&self, d: f64) -> Aabb {
        Aabb {
            min: self.min.map(|x| x - d),
            max: self.max.map(|x| x + d),
        }
    }

    pub fn is_empty(&self) -> bool {
        (0..3).any(|k| self.min[k] > self.max[k])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn box_basics() {
        let b = Aabb::from_points(&[[0.0, 1.0, 2.0], [2.0, -1.0, 3.0]]);
        assert_eq!(
            b,
            Aabb {
                min: [0.0, -1.0, 2.0],
                max: [2.0, 1.0, 3.0]
            }
        );
        assert!(!b.is_empty() && Aabb::empty().is_empty());
        assert!(b.contains([2.0, 1.0, 3.0]) && !b.contains([2.1, 0.0, 2.5]));
        assert_eq!(b.center(), [1.0, 0.0, 2.5]);
        assert_eq!(b.distance2([4.0, 0.0, 2.5]), 4.0);
        assert_eq!(b.distance2([1.0, 0.0, 2.5]), 0.0);
        let touching = Aabb {
            min: [2.0, 1.0, 3.0],
            max: [5.0; 3],
        };
        assert!(b.intersects(&touching));
        assert!(!b.intersects(&touching.expand(-0.5)));
        assert_eq!(Aabb::empty().union(&b), b);
        assert!(!Aabb::empty().intersects(&b));
    }
}
