/// SPHINCS+ address structure (ADRS).
///
/// A 32-byte structure encoding layer, tree address, type, and
/// type-specific fields (keypair, chain, hash, tree height/index).

/// Address types
pub const WOTS_HASH: u32 = 0;
pub const WOTS_PK: u32 = 1;
pub const TREE: u32 = 2;
pub const FORS_TREE: u32 = 3;
pub const FORS_ROOTS: u32 = 4;

/// 32-byte ADRS structure, laid out as 8 x u32 fields:
///   [0]: layer address
///   [1..3]: tree address (64-bit, split across words 1-2; word 3 low)
///   [3]: type
///   [4..7]: type-specific
#[derive(Clone, Debug)]
pub struct Address {
    data: [u8; 32],
}

impl Address {
    pub fn new() -> Self {
        Address { data: [0u8; 32] }
    }

    pub fn as_bytes(&self) -> &[u8; 32] {
        &self.data
    }

    fn set_u32(&mut self, offset: usize, value: u32) {
        self.data[offset..offset + 4].copy_from_slice(&value.to_be_bytes());
    }

    fn get_u32(&self, offset: usize) -> u32 {
        u32::from_be_bytes([
            self.data[offset],
            self.data[offset + 1],
            self.data[offset + 2],
            self.data[offset + 3],
        ])
    }

    pub fn set_layer_address(&mut self, layer: u32) {
        self.set_u32(0, layer);
    }

    pub fn get_layer_address(&self) -> u32 {
        self.get_u32(0)
    }

    pub fn set_tree_address(&mut self, tree: u64) {
        // Tree address occupies bytes 4..12 (words 1 and 2)
        self.data[4..12].copy_from_slice(&tree.to_be_bytes());
    }

    pub fn get_tree_address(&self) -> u64 {
        let mut buf = [0u8; 8];
        buf.copy_from_slice(&self.data[4..12]);
        u64::from_be_bytes(buf)
    }

    pub fn set_type(&mut self, addr_type: u32) {
        self.set_u32(12, addr_type);
        // Clear type-specific fields when type changes
        for i in 16..32 {
            self.data[i] = 0;
        }
    }

    pub fn get_type(&self) -> u32 {
        self.get_u32(12)
    }

    pub fn set_keypair_address(&mut self, keypair: u32) {
        self.set_u32(16, keypair);
    }

    pub fn get_keypair_address(&self) -> u32 {
        self.get_u32(16)
    }

    pub fn set_chain_address(&mut self, chain: u32) {
        self.set_u32(20, chain);
    }

    pub fn set_hash_address(&mut self, hash: u32) {
        self.set_u32(24, hash);
    }

    pub fn set_tree_height(&mut self, height: u32) {
        self.set_u32(20, height);
    }

    pub fn get_tree_height(&self) -> u32 {
        self.get_u32(20)
    }

    pub fn set_tree_index(&mut self, index: u32) {
        self.set_u32(24, index);
    }

    pub fn get_tree_index(&self) -> u32 {
        self.get_u32(24)
    }
}

impl Default for Address {
    fn default() -> Self {
        Self::new()
    }
}
