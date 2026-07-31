<!-- TagCloud.vue -->
<template>
  <div class="card shadow-sm border-0 mb-3" :style="{ backgroundColor: 'var(--app-card-bg)' }">
    <div class="card-header border-0 px-4 pt-3 pb-2" :style="{ backgroundColor: 'var(--app-card-bg)' }">
      <h6 class="m-0" :style="{ color: 'var(--app-text-muted)' }">{{ tagsText }}</h6>
    </div>
    <div class="card-body px-4 pt-2 pb-3">
      <div class="d-flex flex-wrap gap-2">
        <span v-for="t in cloudTags" :key="t.name" class="badge tag-badge fw-normal py-1 px-2 rounded-3 cursor-pointer"
          @click="goTag(t.name)">
          # {{ t.name }}
        </span>
      </div>
    </div>
  </div>
</template>

<script setup>
import { computed, nextTick, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { useRoute, useRouter } from 'vue-router'

const props = defineProps({
  tagData: { type: Array, default: () => [] }
})

const { t } = useI18n()
const route = useRoute()
const router = useRouter()

const tags = ref([])

const tagsText = computed(() => t('tags'))
const cloudTags = computed(() => {
  const src = Array.isArray(props.tagData) && props.tagData.length ? props.tagData : tags.value
  const map = new Map()
  src.forEach(t => {
    const name = typeof t === 'string' ? t : t.name
    if (!name) return
    map.set(name, (map.get(name) || 0) + 1)
  })
  return Array.from(map.entries()).map(([name, count]) => ({ name, count })).sort((a, b) => a.name.localeCompare(b.name))
})

if (props.tagData.length) tags.value = props.tagData

function goTag(name) {
  if (!name) return
  const q = { ...route.query, tag: name, page: '1' }
  router.push({ path: '/', query: q }).catch(() => { })
  nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }))
}
</script>

<style scoped>
.card .badge {
  font-size: 0.85rem;
  font-weight: 500;
}

.tag-badge {
  color: var(--app-tag-text) !important;
  background-color: var(--app-tag-bg) !important;
}
</style>
