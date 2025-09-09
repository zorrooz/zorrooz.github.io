<!-- TagCloud.vue -->
<template>
  <div class="card shadow-sm border-0 bg-white mb-3">
    <div class="card-header bg-white border-0 pt-3 px-4 pb-2">
      <h6 class="text-secondary m-0 small">标签</h6>
    </div>
    <div class="card-body p-4">
      <div class="d-flex flex-wrap gap-2">
        <span v-for="tag in tags" :key="tag.id" class="badge bg-light text-dark fw-normal py-1 px-2 small">
          # {{ tag.name }}
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
      tags: [
        { id: 1, name: 'Vue.js' },
        { id: 2, name: 'JavaScript' },
        { id: 3, name: 'Vue.js' },
        { id: 4, name: 'JavaScript' },
        { id: 5, name: 'Vue.js' },
        { id: 6, name: 'JavaScript' },
        { id: 7, name: 'Vue.js' },
        { id: 8, name: 'JavaScript' },
        { id: 9, name: 'Vue.js' },
        { id: 10, name: 'JavaScript' },
        // 其他标签可扩展
      ]
    }
  },
  mounted() {
    if (this.tagData.length) {
      this.tags = this.tagData
    }
  }
}
</script>

<style scoped></style>
