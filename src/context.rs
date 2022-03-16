use std::fmt;

use crate::depth::BASES;

pub const N_CTXT: usize = 64;

// 5 base sequence context
#[derive(Clone, Copy, Debug)]
pub(crate) struct Ctxt5(pub(crate) u16);

impl fmt::Display for Ctxt5 {
   fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
      let (a, b) = if f.alternate() {
         (1, 4)
      } else {
         (0, 5)
      };
      let x = self.0 as usize;
      for i in a..b {
         let j = (4 - i) << 1;
         write!(f, "{}", BASES.as_bytes()[(x >> j) & 3] as char)?;
      }
      Ok(())
   }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Context5 {
   x: u16,
}

impl fmt::Display for Context5 {
   fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
      if let Some(x) = self.context() {
         let ct = Ctxt5(x as u16);
         if f.alternate() {
            write!(f, "{:#}", ct)?;
         } else {
            write!(f, "{}", ct)?;
         }
      } else {
         write!(f, "{}", if f.alternate() { "---" } else { "-----"} )?;
      }
      Ok(())
   }
}

impl Context5 {
   pub(crate) fn new(reverse: bool) -> Self {
      let x = if reverse { 0x8000 } else { 0 };
      Self { x }
   }

   pub(crate) fn add_base(&mut self, c: u8) {
      let c = (c & 3) as u16;
      if (self.x & 0x8000) == 0 {
         // base context
         let c1 = ((self.x & 0xff) << 2) | c;
         // valid mask
         let c2 = ((self.x & 0x3c00) << 1) | 0x400;
         self.x = c1 | c2;
      } else {
         // base context
         let c1 = ((self.x & 0x3ff) >> 2) | ((c ^ 3) << 8);
         // valid mask
         let c2 = ((self.x & 0x7800) >> 1) | 0x4000;
         self.x = 0x8000 | c1 | c2;
      }
   }

   pub(crate) fn context(&self) -> Option<usize> {
      if (self.x & 0x7c00) == 0x7c00 {
         Some((self.x & 0x3ff) as usize)
      } else {
         None
      }
   }

   pub(crate) fn context3(&self) -> Option<usize> {
      if (self.x & 0x7c00) == 0x7c00 {
         Some(((self.x & 0xff) >> 2) as usize)
      } else {
         None
      }
   }
}

