<!-- TagCloud.vue -->
<template>
  <div class="tag-cloud-card">
    <div class="tag-cloud-header">
      <h6 class="m-0">{{ tagsText }}</h6>
    </div>
    <div class="tag-cloud-body">
      <span v-for="t in cloudTags" :key="t.name" class="tag" @click="goTag(t.name)">
        <i class="fas fa-hashtag"></i>{{ t.name }}
      </span>
    </div>
  </div>
</template>

<script setup lang="ts">
import { computed, nextTick, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { useRoute, useRouter } from 'vue-router'

const props = withDefaults(defineProps<{ tagData?: unknown[] }>(), { tagData: () => [] })

const { t } = useI18n()
const route = useRoute()
const router = useRouter()

const tags = ref<unknown[]>([])

const tagsText = computed(() => t('tags'))
const cloudTags = computed(() => {
  const src = Array.isArray(props.tagData) && props.tagData.length ? props.tagData : tags.value
  const map = new Map<string, number>()
  src.forEach(t => {
    const name = typeof t === 'string' ? t : (t as { name?: string })?.name
    if (!name) return
    map.set(name, (map.get(name) || 0) + 1)
  })
  return Array.from(map.entries()).map(([name, count]) => ({ name, count })).sort((a, b) => a.name.localeCompare(b.name))
})

if (props.tagData.length) tags.value = props.tagData

function goTag(name: string) {
  if (!name) return
  const q = { ...route.query, tag: name, page: '1' }
  router.push({ path: '/', query: q }).catch(() => { })
  nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }))
}
</script>

<style scoped>
.tag-cloud-card {
  background-color: var(--surface);
  border: 1px solid var(--line);
  border-radius: var(--radius);
  transition: border-color 0.22s ease, box-shadow 0.22s ease;
}

.tag-cloud-card:hover {
  border-color: color-mix(in srgb, var(--primary) 30%, transparent);
  box-shadow: var(--shadow-soft);
}

.tag-cloud-header {
  padding: 1.25rem 1.5rem 0.5rem;
}

.tag-cloud-header h6 {
  font-size: var(--text-base);
  font-weight: 700;
  letter-spacing: -0.01em;
  color: var(--fg);
}

.tag-cloud-body {
  display: flex;
  flex-wrap: wrap;
  gap: var(--sp-2);
  padding: 0.5rem 1.5rem 1.25rem;
}
</style>
