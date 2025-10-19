<!-- TagCloud.vue -->
<template>
  <div class="card shadow-sm border-0 mb-3" :style="{ backgroundColor: 'var(--app-card-bg)' }">
    <div class="card-header border-0 px-4 py-2" :style="{ backgroundColor: 'var(--app-card-bg)' }">
      <h6 class="m-0" :style="{ color: 'var(--app-text-muted)' }">标签</h6>
    </div>
    <div class="card-body px-4 pt-2 pb-3">
      <div class="d-flex flex-wrap gap-2">
        <span v-for="t in cloudTags" :key="t.name" class="badge tag-badge fw-normal py-1 px-2 rounded-3" style="cursor: pointer" @click="goTag(t.name)">
          # {{ t.name }}
        </span>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  name: 'TagCloud',
  props: {
    tagData: {
      type: Array,
      default: () => []
    }
  },
  data() {
    return {
      tags: []
    }
  },
  computed: {
    cloudTags() {
      const src = Array.isArray(this.tagData) && this.tagData.length ? this.tagData : this.tags;
      const map = new Map();
      src.forEach(t => {
        const name = typeof t === 'string' ? t : t.name;
        if (!name) return;
        map.set(name, (map.get(name) || 0) + 1);
      });
      return Array.from(map.entries()).map(([name, count]) => ({ name, count }))
        .sort((a, b) => a.name.localeCompare(b.name));
    }
  },
  mounted() {
    if (this.tagData.length) {
      this.tags = this.tagData
    }
  },
  methods: {
    goTag(name) {
      if (!name) return;
      const q = { ...this.$route.query, tag: name, page: '1' };
      this.$router.push({ path: '/', query: q }).catch(() => {});
      this.$nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }));
    }
  }
}
</script>

<style scoped>
.card .badge { font-size: 0.85rem; font-weight: 500; }
.tag-badge { 
  color: var(--app-tag-text) !important; 
  background-color: var(--app-tag-bg) !important; 
}
</style>
